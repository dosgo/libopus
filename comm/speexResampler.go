package comm

import (
	"math"
	"sync"
)

const (
	fixedStackAlloc = 8192
)

// SpeexResampler 实现任意比率的音频重采样
type SpeexResampler struct {
	mu sync.Mutex

	inRate        int
	outRate       int
	numRate       int
	denRate       int
	quality       int
	nbChannels    int
	filtLen       int
	memAllocSize  int
	bufferSize    int
	intAdvance    int
	fracAdvance   int
	cutoff        float64
	oversample    int
	initialised   int
	started       int
	lastSample    []int
	sampFracNum   []int
	magicSamples  []int
	mem           []float64
	sincTable     []float64
	sincTableLen  int
	resamplerFunc func(channelIndex int, input []float64, inLen int, output []float64, outLen int) (int, int)
	inStride      int
	outStride     int
}

// NewSpeexResampler 创建新的重采样器
func NewSpeexResampler(nbChannels, inRate, outRate, quality int) *SpeexResampler {
	return NewSpeexResamplerFractional(nbChannels, inRate, outRate, inRate, outRate, quality)
}

// NewSpeexResamplerFractional 创建支持分数比率的重采样器
func NewSpeexResamplerFractional(nbChannels, ratioNum, ratioDen, inRate, outRate, quality int) *SpeexResampler {
	if quality > 10 || quality < 0 {
		panic("quality must be between 0 and 10")
	}

	r := &SpeexResampler{
		nbChannels:   nbChannels,
		quality:      -1,
		bufferSize:   160,
		inStride:     1,
		outStride:    1,
		lastSample:   make([]int, nbChannels),
		magicSamples: make([]int, nbChannels),
		sampFracNum:  make([]int, nbChannels),
	}

	r.SetQuality(quality)
	r.SetRateFraction(ratioNum, ratioDen, inRate, outRate)
	r.updateFilter()
	r.initialised = 1
	return r
}

// 质量映射表
type qualityMapping struct {
	baseLength          int
	oversample          int
	downsampleBandwidth float64
	upsampleBandwidth   float64
	windowFunc          *funcDef
}

type funcDef struct {
	table      []float64
	oversample int
}

var qualityMap = []qualityMapping{
	{8, 4, 0.830, 0.860, &funcDef{kaiser6Table, 32}},
	{16, 4, 0.850, 0.880, &funcDef{kaiser6Table, 32}},
	{32, 4, 0.882, 0.910, &funcDef{kaiser6Table, 32}},
	{48, 8, 0.895, 0.917, &funcDef{kaiser8Table, 32}},
	{64, 8, 0.921, 0.940, &funcDef{kaiser8Table, 32}},
	{80, 16, 0.922, 0.940, &funcDef{kaiser10Table, 32}},
	{96, 16, 0.940, 0.945, &funcDef{kaiser10Table, 32}},
	{128, 16, 0.950, 0.950, &funcDef{kaiser10Table, 32}},
	{160, 16, 0.960, 0.960, &funcDef{kaiser10Table, 32}},
	{192, 32, 0.968, 0.968, &funcDef{kaiser12Table, 64}},
	{256, 32, 0.975, 0.975, &funcDef{kaiser12Table, 64}},
}

var (
	kaiser12Table = []float64{
		0.99859849, 1.00000000, 0.99859849, 0.99440475, 0.98745105, 0.97779076,
		0.96549770, 0.95066529, 0.93340547, 0.91384741, 0.89213598, 0.86843014,
		0.84290116, 0.81573067, 0.78710866, 0.75723148, 0.72629970, 0.69451601,
		0.66208321, 0.62920216, 0.59606986, 0.56287762, 0.52980938, 0.49704014,
		0.46473455, 0.43304576, 0.40211431, 0.37206735, 0.34301800, 0.31506490,
		0.28829195, 0.26276832, 0.23854851, 0.21567274, 0.19416736, 0.17404546,
		0.15530766, 0.13794294, 0.12192957, 0.10723616, 0.09382272, 0.08164178,
		0.07063950, 0.06075685, 0.05193064, 0.04409466, 0.03718069, 0.03111947,
		0.02584161, 0.02127838, 0.01736250, 0.01402878, 0.01121463, 0.00886058,
		0.00691064, 0.00531256, 0.00401805, 0.00298291, 0.00216702, 0.00153438,
		0.00105297, 0.00069463, 0.00043489, 0.00025272, 0.00013031, 0.0000527734,
		0.00001000, 0.00000000,
	}

	kaiser10Table = []float64{
		0.99537781, 1.00000000, 0.99537781, 0.98162644, 0.95908712, 0.92831446,
		0.89005583, 0.84522401, 0.79486424, 0.74011713, 0.68217934, 0.62226347,
		0.56155915, 0.50119680, 0.44221549, 0.38553619, 0.33194107, 0.28205962,
		0.23636152, 0.19515633, 0.15859932, 0.12670280, 0.09935205, 0.07632451,
		0.05731132, 0.04193980, 0.02979584, 0.02044510, 0.01345224, 0.00839739,
		0.00488951, 0.00257636, 0.00115101, 0.00035515, 0.00000000, 0.00000000,
	}

	kaiser8Table = []float64{
		0.99635258, 1.00000000, 0.99635258, 0.98548012, 0.96759014, 0.94302200,
		0.91223751, 0.87580811, 0.83439927, 0.78875245, 0.73966538, 0.68797126,
		0.63451750, 0.58014482, 0.52566725, 0.47185369, 0.41941150, 0.36897272,
		0.32108304, 0.27619388, 0.23465776, 0.19672670, 0.16255380, 0.13219758,
		0.10562887, 0.08273982, 0.06335451, 0.04724088, 0.03412321, 0.02369490,
		0.01563093, 0.00959968, 0.00527363, 0.00233883, 0.00050000, 0.00000000,
	}

	kaiser6Table = []float64{
		0.99733006, 1.00000000, 0.99733006, 0.98935595, 0.97618418, 0.95799003,
		0.93501423, 0.90755855, 0.87598009, 0.84068475, 0.80211977, 0.76076565,
		0.71712752, 0.67172623, 0.62508937, 0.57774224, 0.53019925, 0.48295561,
		0.43647969, 0.39120616, 0.34752997, 0.30580127, 0.26632152, 0.22934058,
		0.19505503, 0.16360756, 0.13508755, 0.10953262, 0.08693120, 0.06722600,
		0.05031820, 0.03607231, 0.02432151, 0.01487334, 0.00752000, 0.00000000,
	}
)

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func (r *SpeexResampler) float2int(x float64) int16 {
	if x < -32768 {
		return -32768
	}
	if x > 32767 {
		return 32767
	}
	return int16(x)
}

func computeFunc(x float64, f *funcDef) float64 {
	y := float64(x) * float64(f.oversample)
	ind := int(math.Floor(y))
	frac := y - float64(ind)

	interp3 := -0.1666666667*frac + 0.1666666667*math.Pow(frac, 3)
	interp2 := frac + 0.5*math.Pow(frac, 2) - 0.5*math.Pow(frac, 3)
	interp0 := -0.3333333333*frac + 0.5*math.Pow(frac, 2) - 0.1666666667*math.Pow(frac, 3)
	interp1 := 1.0 - interp3 - interp2 - interp0

	return interp0*f.table[ind] + interp1*f.table[ind+1] + interp2*f.table[ind+2] + interp3*f.table[ind+3]
}

func sinc(cutoff, x float64, N int, window *funcDef) float64 {
	/*fprintf (stderr, "%f ", x);*/
	xx := x * cutoff
	if math.Abs(x) < 1e-6 {
		return cutoff
	} else if math.Abs(x) > 0.5*float64(N) {
		return 0
	}
	/*FIXME: Can it really be any slower than this? */
	return float64(cutoff * math.Sin(math.Pi*xx) / (math.Pi * xx) * computeFunc(math.Abs(2.0*x/float64(N)), window))

}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func cubicCoef(frac float64, interp []float64) {
	interp[0] = -0.16667*frac + 0.16667*frac*frac*frac
	interp[1] = frac + 0.5*frac*frac - 0.5*frac*frac*frac
	interp[3] = -0.33333*frac + 0.5*frac*frac - 0.16667*frac*frac*frac
	interp[2] = 1.0 - interp[0] - interp[1] - interp[3]
}

func (r *SpeexResampler) resamplerBasicDirectSingle(channelIndex int, input []float64, inLen int, output []float64, outLen int) (int, int) {
	N := r.filtLen
	outSample := 0
	lastSample := r.lastSample[channelIndex]
	sampFracNum := r.sampFracNum[channelIndex]

	for lastSample < inLen && outSample < outLen {
		sinct := sampFracNum * N
		sum := float64(0)
		for j := 0; j < N; j++ {
			sum += r.sincTable[sinct+j] * input[lastSample+j]
		}

		output[outSample] = sum
		outSample++
		lastSample += r.intAdvance
		sampFracNum += r.fracAdvance
		if sampFracNum >= r.denRate {
			sampFracNum -= r.denRate
			lastSample++
		}
	}

	r.lastSample[channelIndex] = lastSample
	r.sampFracNum[channelIndex] = sampFracNum
	return outSample, lastSample
}

func (r *SpeexResampler) resamplerBasicInterpolateSingle(channelIndex int, input []float64, inLen int, output []float64, outLen int) (int, int) {
	N := r.filtLen
	outSample := 0
	lastSample := r.lastSample[channelIndex]
	sampFracNum := r.sampFracNum[channelIndex]

	interp := make([]float64, 4)
	accum := make([]float64, 4)

	for lastSample < inLen && outSample < outLen {
		offset := sampFracNum * r.oversample / r.denRate
		frac := float64((sampFracNum*r.oversample)%r.denRate) / float64(r.denRate)

		for i := range accum {
			accum[i] = 0
		}

		for j := 0; j < N; j++ {
			currIn := input[lastSample+j]
			accum[0] += currIn * r.sincTable[4+(j+1)*r.oversample-offset-2]
			accum[1] += currIn * r.sincTable[4+(j+1)*r.oversample-offset-1]
			accum[2] += currIn * r.sincTable[4+(j+1)*r.oversample-offset]
			accum[3] += currIn * r.sincTable[4+(j+1)*r.oversample-offset+1]
		}

		cubicCoef(frac, interp)
		sum := interp[0]*accum[0] + interp[1]*accum[1] + interp[2]*accum[2] + interp[3]*accum[3]

		output[outSample] = sum
		outSample++
		lastSample += r.intAdvance
		sampFracNum += r.fracAdvance
		if sampFracNum >= r.denRate {
			sampFracNum -= r.denRate
			lastSample++
		}
	}

	r.lastSample[channelIndex] = lastSample
	r.sampFracNum[channelIndex] = sampFracNum
	return outSample, lastSample
}

func (r *SpeexResampler) updateFilter() {
	r.mu.Lock()
	defer r.mu.Unlock()

	oldLength := r.filtLen
	qmap := qualityMap[r.quality]
	r.oversample = qmap.oversample
	r.filtLen = qmap.baseLength

	if r.numRate > r.denRate {
		r.cutoff = qmap.downsampleBandwidth * float64(r.denRate) / float64(r.numRate)
		r.filtLen = r.filtLen * r.numRate / r.denRate
		r.filtLen = ((r.filtLen - 1) &^ 0x7) + 8
		if 2*r.denRate < r.numRate {
			r.oversample >>= 1
		}
		if 4*r.denRate < r.numRate {
			r.oversample >>= 1
		}
		if 8*r.denRate < r.numRate {
			r.oversample >>= 1
		}
		if 16*r.denRate < r.numRate {
			r.oversample >>= 1
		}
		if r.oversample < 1 {
			r.oversample = 1
		}
	} else {
		r.cutoff = qmap.upsampleBandwidth
	}

	if r.denRate <= 16*(r.oversample+8) {
		tableSize := r.filtLen * r.denRate
		if len(r.sincTable) < tableSize {
			r.sincTable = make([]float64, tableSize)
		}

		for i := 0; i < r.denRate; i++ {
			for j := 0; j < r.filtLen; j++ {
				pos := j - r.filtLen/2 + 1
				r.sincTable[i*r.filtLen+j] = sinc(
					r.cutoff,
					float64(pos)-float64(i)/float64(r.denRate),
					r.filtLen,
					qmap.windowFunc,
				)
			}
		}
		r.resamplerFunc = r.resamplerBasicDirectSingle
	} else {
		tableSize := r.filtLen*r.oversample + 8
		if len(r.sincTable) < tableSize {
			r.sincTable = make([]float64, tableSize)
		}

		for i := -4; i < r.oversample*r.filtLen+4; i++ {
			idx := i + 4
			pos := float64(i)/float64(r.oversample) - float64(r.filtLen)/2
			r.sincTable[idx] = sinc(r.cutoff, pos, r.filtLen, qmap.windowFunc)
		}
		r.resamplerFunc = r.resamplerBasicInterpolateSingle
	}

	r.intAdvance = r.numRate / r.denRate
	r.fracAdvance = r.numRate % r.denRate

	if r.mem == nil {
		r.memAllocSize = r.filtLen - 1 + r.bufferSize
		r.mem = make([]float64, r.nbChannels*r.memAllocSize)
	} else if r.started == 0 {
		r.memAllocSize = r.filtLen - 1 + r.bufferSize
		r.mem = make([]float64, r.nbChannels*r.memAllocSize)
	} else if r.filtLen > oldLength {
		oldAllocSize := r.memAllocSize
		if r.filtLen-1+r.bufferSize > r.memAllocSize {
			r.memAllocSize = r.filtLen - 1 + r.bufferSize
			r.mem = make([]float64, r.nbChannels*r.memAllocSize)
		}

		for i := 0; i < r.nbChannels; i++ {
			offset := i * r.memAllocSize
			oldOffset := i * oldAllocSize
			olen := oldLength

			if r.magicSamples[i] > 0 {
				olen = oldLength + 2*r.magicSamples[i]
				for j := oldLength - 2 + r.magicSamples[i]; j >= 0; j-- {
					r.mem[offset+j+r.magicSamples[i]] = r.mem[oldOffset+j]
				}
				for j := 0; j < r.magicSamples[i]; j++ {
					r.mem[offset+j] = 0
				}
				r.magicSamples[i] = 0
			}

			if r.filtLen > olen {
				for j := 0; j < olen-1; j++ {
					r.mem[offset+(r.filtLen-2-j)] = r.mem[offset+(olen-2-j)]
				}
				for j := olen - 1; j < r.filtLen-1; j++ {
					r.mem[offset+(r.filtLen-2-j)] = 0
				}
				r.lastSample[i] += (r.filtLen - olen) / 2
			} else {
				r.magicSamples[i] = (olen - r.filtLen) / 2
				for j := 0; j < r.filtLen-1+r.magicSamples[i]; j++ {
					r.mem[offset+j] = r.mem[offset+j+r.magicSamples[i]]
				}
			}
		}
	} else if r.filtLen < oldLength {
		for i := 0; i < r.nbChannels; i++ {
			oldMagic := r.magicSamples[i]
			r.magicSamples[i] = (oldLength - r.filtLen) / 2
			offset := i * r.memAllocSize

			for j := 0; j < r.filtLen-1+r.magicSamples[i]+oldMagic; j++ {
				r.mem[offset+j] = r.mem[offset+j+r.magicSamples[i]]
			}
			r.magicSamples[i] += oldMagic
		}
	}
}

func (r *SpeexResampler) speexResamplerProcessNative(channelIndex int, inLen int, output []float64, outLen int) (int, int) {
	N := r.filtLen
	memPtr := channelIndex * r.memAllocSize
	ilen := inLen

	outSample, lastSample := r.resamplerFunc(channelIndex, r.mem[memPtr:], ilen, output, outLen)

	if lastSample < inLen {
		ilen = lastSample
	}
	outLen = outSample
	r.lastSample[channelIndex] -= ilen

	for j := memPtr; j < N-1+memPtr; j++ {
		r.mem[j] = r.mem[j+ilen]
	}

	return ilen, outLen
}

func (r *SpeexResampler) speexResamplerMagic(channelIndex int, output []float64, outLen int) (int, int) {
	tmpInLen := r.magicSamples[channelIndex]
	outLen, _ = r.speexResamplerProcessNative(channelIndex, tmpInLen, output, outLen)
	r.magicSamples[channelIndex] -= tmpInLen

	if r.magicSamples[channelIndex] != 0 {
		N := r.filtLen
		memPtr := channelIndex * r.memAllocSize
		for j := memPtr; j < r.magicSamples[channelIndex]+memPtr; j++ {
			r.mem[N-1+j] = r.mem[N-1+j+tmpInLen]
		}
	}

	return outLen, outLen
}

// 公共API
func (r *SpeexResampler) ProcessFloat(channelIndex int, input []float64, inLen int, output []float64, outLen int) (int, int) {
	r.mu.Lock()
	defer r.mu.Unlock()

	if channelIndex < 0 || channelIndex >= r.nbChannels {
		panic("invalid channel index")
	}

	ilen := inLen
	olen := outLen
	memPtr := channelIndex * r.memAllocSize
	filtOffs := r.filtLen - 1
	xlen := r.memAllocSize - filtOffs

	if r.magicSamples[channelIndex] != 0 {
		olen, _ = r.speexResamplerMagic(channelIndex, output, olen)
	}

	if r.magicSamples[channelIndex] == 0 {
		for ilen > 0 && olen > 0 {
			ichunk := minInt(ilen, xlen)
			ochunk := olen
			var outProcessed, inProcessed int

			if input != nil {
				copy(r.mem[memPtr+filtOffs:], input[:ichunk])
			} else {
				for j := 0; j < ichunk; j++ {
					r.mem[memPtr+filtOffs+j] = 0
				}
			}

			inProcessed, outProcessed = r.speexResamplerProcessNative(channelIndex, ichunk, output, ochunk)
			ilen -= inProcessed
			olen -= outProcessed
			output = output[outProcessed:]
			if input != nil {
				input = input[ichunk:]
			}
		}
	}

	return inLen - ilen, outLen - olen
}

func (r *SpeexResampler) ProcessShort(channelIndex int, input []int16, inLen int, output []int16, outLen int) (int, int) {
	r.mu.Lock()
	defer r.mu.Unlock()

	if channelIndex < 0 || channelIndex >= r.nbChannels {
		panic("invalid channel index")
	}

	ilen := inLen
	olen := outLen
	memPtr := channelIndex * r.memAllocSize
	xlen := r.memAllocSize - (r.filtLen - 1)

	ystack := make([]float64, minInt(olen, fixedStackAlloc))

	for ilen > 0 && olen > 0 {
		ichunk := minInt(ilen, xlen)
		ochunk := minInt(olen, len(ystack))
		omagic := 0
		y := 0

		if r.magicSamples[channelIndex] != 0 {
			processed, _ := r.speexResamplerMagic(channelIndex, ystack, ochunk)
			omagic = processed
			ochunk -= omagic
			olen -= omagic
		}

		if r.magicSamples[channelIndex] == 0 {
			if input != nil {
				for j := 0; j < ichunk; j++ {
					r.mem[memPtr+j+r.filtLen-1] = float64(input[j])
				}
			} else {
				for j := 0; j < ichunk; j++ {
					r.mem[memPtr+j+r.filtLen-1] = 0
				}
			}

			_, outProcessed := r.speexResamplerProcessNative(channelIndex, ichunk, ystack[y:], ochunk)
			ichunk = outProcessed
			ochunk = outProcessed
		}

		for j := 0; j < ochunk+omagic; j++ {
			output[j] = r.float2int(ystack[j])
		}

		ilen -= ichunk
		olen -= ochunk
		output = output[ochunk+omagic:]
		if input != nil {
			input = input[ichunk:]
		}
	}

	return inLen - ilen, outLen - olen
}

func (r *SpeexResampler) ProcessInterleavedFloat(input []float64, inLen int, output []float64, outLen int) (int, int) {
	origInLen := inLen
	origOutLen := outLen
	inStride := r.inStride
	outStride := r.outStride
	r.inStride = r.nbChannels
	r.outStride = r.nbChannels

	for i := 0; i < r.nbChannels; i++ {
		var chanIn, chanOut []float64
		if input != nil {
			chanIn = input[i:]
		}
		chanOut = output[i:]
		processedIn, processedOut := r.ProcessFloat(i, chanIn, inLen, chanOut, outLen)
		if processedIn < inLen {
			inLen = processedIn
		}
		if processedOut < outLen {
			outLen = processedOut
		}
	}

	r.inStride = inStride
	r.outStride = outStride
	return origInLen - inLen, origOutLen - outLen
}

func (r *SpeexResampler) ProcessInterleavedShort(input []int16, inLen int, output []int16, outLen int) (int, int) {
	origInLen := inLen
	origOutLen := outLen
	inStride := r.inStride
	outStride := r.outStride
	r.inStride = r.nbChannels
	r.outStride = r.nbChannels

	for i := 0; i < r.nbChannels; i++ {
		var chanIn, chanOut []int16
		if input != nil {
			chanIn = input[i:]
		}
		chanOut = output[i:]
		processedIn, processedOut := r.ProcessShort(i, chanIn, inLen, chanOut, outLen)
		if processedIn < inLen {
			inLen = processedIn
		}
		if processedOut < outLen {
			outLen = processedOut
		}
	}

	r.inStride = inStride
	r.outStride = outStride
	return origInLen - inLen, origOutLen - outLen
}

func (r *SpeexResampler) SkipZeroes() {
	for i := 0; i < r.nbChannels; i++ {
		r.lastSample[i] = r.filtLen / 2
	}
}

func (r *SpeexResampler) ResetMem() {
	for i := 0; i < r.nbChannels; i++ {
		r.lastSample[i] = 0
		r.magicSamples[i] = 0
		r.sampFracNum[i] = 0
	}
	for i := range r.mem {
		r.mem[i] = 0
	}
}

// Getters and Setters
func (r *SpeexResampler) SetRates(inRate, outRate int) {
	r.SetRateFraction(inRate, outRate, inRate, outRate)
}

func (r *SpeexResampler) GetRates() (int, int) {
	return r.inRate, r.outRate
}

func (r *SpeexResampler) SetRateFraction(ratioNum, ratioDen, inRate, outRate int) {
	if r.inRate == inRate && r.outRate == outRate && r.numRate == ratioNum && r.denRate == ratioDen {
		return
	}

	oldDen := r.denRate
	r.inRate = inRate
	r.outRate = outRate
	r.numRate = ratioNum
	r.denRate = ratioDen

	// 简化分数
	gcd := func(a, b int) int {
		for b != 0 {
			a, b = b, a%b
		}
		return a
	}

	div := gcd(r.numRate, r.denRate)
	r.numRate /= div
	r.denRate /= div

	if oldDen > 0 {
		for i := 0; i < r.nbChannels; i++ {
			r.sampFracNum[i] = r.sampFracNum[i] * r.denRate / oldDen
			if r.sampFracNum[i] >= r.denRate {
				r.sampFracNum[i] = r.denRate - 1
			}
		}
	}

	if r.initialised != 0 {
		r.updateFilter()
	}
}

func (r *SpeexResampler) GetRateFraction() (int, int) {
	return r.numRate, r.denRate
}

func (r *SpeexResampler) SetQuality(quality int) {
	if quality > 10 || quality < 0 {
		panic("quality must be between 0 and 10")
	}
	if r.quality == quality {
		return
	}
	r.quality = quality
	if r.initialised != 0 {
		r.updateFilter()
	}
}

func (r *SpeexResampler) GetQuality() int {
	return r.quality
}

func (r *SpeexResampler) SetInputStride(stride int) {
	r.inStride = stride
}

func (r *SpeexResampler) SetOutputStride(stride int) {
	r.outStride = stride
}

func (r *SpeexResampler) InputLatency() int {
	return r.filtLen / 2
}

func (r *SpeexResampler) OutputLatency() int {
	return ((r.filtLen/2)*r.denRate + (r.numRate >> 1)) / r.numRate
}
