package org.concentus;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import org.concentus.*;


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author lostromb
 */
public class Program {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
       test5();
    }
    
    public static void test()
    {
        try {
            FileInputStream fileIn = new FileInputStream("./48Khz Stereo.raw");
            OpusEncoder encoder = new OpusEncoder(48000, 2, OpusApplication.OPUS_APPLICATION_AUDIO);
            encoder.setBitrate(96000);
            encoder.setForceMode(OpusMode.MODE_SILK_ONLY);
            encoder.setSignalType(OpusSignal.OPUS_SIGNAL_MUSIC);
            encoder.setComplexity(0);
            
            OpusDecoder decoder = new OpusDecoder(48000, 2);

            FileOutputStream fileOut = new FileOutputStream("./out_j.raw");
            int packetSamples = 960;
            byte[] inBuf = new byte[packetSamples * 2 * 2];
            byte[] data_packet = new byte[1275];
            long start = System.currentTimeMillis();
            int i=0;
            while (fileIn.available() >= inBuf.length) {
            
                int bytesRead = fileIn.read(inBuf, 0, inBuf.length);
                short[] pcm = BytesToShorts(inBuf, 0, inBuf.length);
           
                if(i>1000){
                 break;
                }
                       
              
                     System.out.println("imput md5:" +Arrays.generateMD5(inBuf));
                int bytesEncoded = encoder.encode(pcm, 0, packetSamples, data_packet, 0, 1275);
                if(i==4){
                    Arrays.debug = false;
                }else{
                    Arrays.debug = false;
                 
                }
             
               
                
                      
                    System.out.println("bytesEncoded:"+bytesEncoded+" data_packet:" +Arrays.generateMD5(data_packet));
        
                    int samplesDecoded = decoder.decode(data_packet, 0, bytesEncoded, pcm, 0, packetSamples, false);
                
                 
                    byte[] bytesOut = ShortsToBytes(pcm);
                               System.out.println("pcm:" + Arrays.generateMD5(pcm)); 
             
                       i++;
           
                   
                
            }
            
            long end = System.currentTimeMillis();
            System.out.println("Time was " + (end - start) + "ms");
            fileIn.close();
            fileOut.close();
            System.out.println("Done!");
        } catch (IOException e) {
            System.out.println(e.getMessage());
        } catch (OpusException e) {
            System.out.println(e.getMessage());
        }
    }


    public static void test1()
    {
        try {
            FileInputStream fileIn = new FileInputStream("./16Khz Mono.raw");
            OpusEncoder encoder = new OpusEncoder(16000, 1, OpusApplication.OPUS_APPLICATION_AUDIO);

            encoder.setBitrate(96000);
            encoder.setForceMode(OpusMode.MODE_CELT_ONLY);
            encoder.setSignalType(OpusSignal.OPUS_SIGNAL_MUSIC);
            encoder.setComplexity(0);

            
            
            OpusDecoder decoder = new OpusDecoder(16000, 1);

            FileOutputStream fileOut = new FileOutputStream("./out_j.raw");
            int packetSamples = 960;
            byte[] inBuf = new byte[packetSamples * 2];
            byte[] data_packet = new byte[1275];
            long start = System.currentTimeMillis();
            int i=0;
            while (fileIn.available() >= inBuf.length) {
            
                int bytesRead = fileIn.read(inBuf, 0, inBuf.length);
                short[] pcm = BytesToShorts(inBuf, 0, inBuf.length);
            
              
                if(i>1000){
                 break;
                }
                         if(i==20){
                        Arrays.debug = false;
                    }else{
                        Arrays.debug = false;
                    }
             
                  
                     System.out.println("imput md5:" +Arrays.generateMD5(inBuf));
                int bytesEncoded = encoder.encode(pcm, 0, packetSamples, data_packet, 0, 1275);
                 
             
               
                
                      
                    System.out.println("bytesEncoded:"+bytesEncoded+" data_packet:" +Arrays.generateMD5(data_packet));
                    int samplesDecoded = decoder.decode(data_packet, 0, bytesEncoded, pcm, 0, packetSamples, false);
        
                 
                    byte[] bytesOut = ShortsToBytes(pcm);
                    System.out.println("pcm:" + Arrays.generateMD5(pcm)); 
                    i++;                
            }
            
            long end = System.currentTimeMillis();
            System.out.println("Time was " + (end - start) + "ms");
            fileIn.close();
            fileOut.close();
            System.out.println("Done!");
        } catch (IOException e) {
            System.out.println(e.getMessage());
        } catch (OpusException e) {
            System.out.println(e.getMessage());
        }
    }


    /// <summary>
    /// Converts interleaved byte samples (such as what you get from a capture device)
    /// into linear short samples (that are much easier to work with)
    /// </summary>
    /// <param name="input"></param>
    /// <returns></returns>
    public static short[] BytesToShorts(byte[] input) {
        return BytesToShorts(input, 0, input.length);
    }

    /// <summary>
    /// Converts interleaved byte samples (such as what you get from a capture device)
    /// into linear short samples (that are much easier to work with)
    /// </summary>
    /// <param name="input"></param>
    /// <returns></returns>
    public static short[] BytesToShorts(byte[] input, int offset, int length) {
        short[] processedValues = new short[length / 2];
        for (int c = 0; c < processedValues.length; c++) {
            short a = (short) (((int) input[(c * 2) + offset]) & 0xFF);
            short b = (short) (((int) input[(c * 2) + 1 + offset]) << 8);
            processedValues[c] = (short) (a | b);
        }

        return processedValues;
    }

    /// <summary>
    /// Converts linear short samples into interleaved byte samples, for writing to a file, waveout device, etc.
    /// </summary>
    /// <param name="input"></param>
    /// <returns></returns>
    public static byte[] ShortsToBytes(short[] input) {
        return ShortsToBytes(input, 0, input.length);
    }

    /// <summary>
    /// Converts linear short samples into interleaved byte samples, for writing to a file, waveout device, etc.
    /// </summary>
    /// <param name="input"></param>
    /// <returns></returns>
    public static byte[] ShortsToBytes(short[] input, int offset, int length) {
        byte[] processedValues = new byte[length * 2];
        for (int c = 0; c < length; c++) {
            processedValues[c * 2] = (byte) (input[c + offset] & 0xFF);
            processedValues[c * 2 + 1] = (byte) ((input[c + offset] >> 8) & 0xFF);
        }

        return processedValues;
    }

        public static void test5()
    {
       try {
            FileInputStream fileIn = new FileInputStream("./testnew.opus");
          OpusDecoder decoder = new OpusDecoder(48000, 2);


            byte[] inBuf = new byte[960*20];
            int i=0;
                while (true) {
            


                 byte[] lengthBytes = new byte[4];
                int bytesReadLen = fileIn.read(lengthBytes);
                
                if (bytesReadLen != 4) {
                   return ;
                }
                
                // 2. 将长度前缀转换为整数（大端序）
                int dataLength = ByteBuffer.wrap(lengthBytes)
                        .order(ByteOrder.LITTLE_ENDIAN)
                        .getInt();
                
                
    

                int bytesRead = fileIn.read(inBuf, 0, dataLength);
                   byte[] partial = java.util.Arrays.copyOfRange(inBuf, 0, Math.min(bytesRead, inBuf.length));
                System.out.printf("inBuf:%s\r\n", java.util.Arrays.toString(partial));
                
                short[] pcm = new short[960*2];
                if (i == 4) {
			        Arrays.debug = true;
		        }

                int samplesDecoded = decoder.decode(inBuf, 0, bytesRead, pcm, 0, 9600, false);
                
            //     System.out.println("pcm:"+java.util.Arrays.toString(pcm));

                  short[] pcmout=new short[samplesDecoded];
                  System.arraycopy(pcm, 0, pcmout, 0, samplesDecoded);

                System.out.printf("pcm:%s\r\n", java.util.Arrays.toString(pcmout));

                
                System.out.println("i:"+i+"len:"+samplesDecoded+"pcm:" + Arrays.generateMD5(pcm,samplesDecoded)); 
                if(fileIn.available() < inBuf.length){
                    break;
                }
                i++;
                    
                if(i>4){
                    break;
                }
            }
            
            

            fileIn.close();
         
            System.out.println("Done!");
             } catch (IOException e) {
            System.out.println(e.getMessage());
        } catch (OpusException e) {
            System.out.println(e.getMessage());
        }
    
    }
}
