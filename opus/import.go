package opus

import (
	"github.com/dosgo/libopus/celt"
	"github.com/dosgo/libopus/comm"
	"github.com/dosgo/libopus/silk"
)

var inlines = comm.Inlines{}
var kernels = comm.Kernels{}
var SilkConstants = silk.SilkConstants
var SilkTables = silk.SilkTables

var TuningParameters = silk.TuningParameters

var SilkError = silk.SilkError

var OpusFramesizeHelpers = celt.OpusFramesizeHelpers
