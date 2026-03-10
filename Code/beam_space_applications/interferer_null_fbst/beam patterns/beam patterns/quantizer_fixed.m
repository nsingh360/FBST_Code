function [A_quant] = quantizer_fixed(A,bits,scale,offset)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% scale and split transmit signal into real and imaginary channels
[M,N] = size(A);

A_real_scaled = (real(A)-offset)/scale;
A_imag_scaled = (imag(A)-offset)/scale;

% quantize DAC signal
DAC_bits = bits;
DAC_bit_step = 2^(-DAC_bits+1);
DAC_codebook = [(-2^(bits-1)):(2^(bits-1)-1)]/2^(bits-1);
DAC_partition = [(-2^(bits-1)+1):(2^(bits-1)-1)]/2^(bits-1);

% quantize real part of DAC
[A_real_index,A_real_quant] = quantiz(A_real_scaled(:),DAC_partition,DAC_codebook);

% quantize imag part of DAC
[A_imag_index,A_imag_quant] = quantiz(A_imag_scaled(:),DAC_partition,DAC_codebook);

A_quant = (reshape(A_real_quant,M,N) + 1i*reshape(A_imag_quant,M,N))*scale + offset + 1i*offset;
end