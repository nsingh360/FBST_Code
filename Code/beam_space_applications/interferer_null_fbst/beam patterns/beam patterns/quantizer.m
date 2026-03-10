function [A_quant,offset_DAC,scale_DAC] = quantizer(A,bits)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% scale and split transmit signal into real and imaginary channels
[M,N] = size(A);
max_DAC = max(max(real(A(:))),max(imag(A(:))));
min_DAC = min(min(real(A(:))),min(imag(A(:))));
scale_DAC = (max_DAC-min_DAC)/2;
offset_DAC = (max_DAC+min_DAC)/2;

A_real_scaled = (real(A)-offset_DAC)/scale_DAC;
A_imag_scaled = (imag(A)-offset_DAC)/scale_DAC;

% quantize DAC signal
DAC_bits = bits;
DAC_bit_step = 2^(-DAC_bits+1);
DAC_codebook = [(-2^(bits-1)):(2^(bits-1)-1)]/2^(bits-1);
DAC_partition = [(-2^(bits-1)+1):(2^(bits-1)-1)]/2^(bits-1);

% quantize real part of DAC
[A_real_index,A_real_quant] = quantiz(A_real_scaled(:),DAC_partition,DAC_codebook);

% quantize imag part of DAC
[A_imag_index,A_imag_quant] = quantiz(A_imag_scaled(:),DAC_partition,DAC_codebook);

A_quant = (reshape(A_real_quant,M,N) + 1i*reshape(A_imag_quant,M,N))*scale_DAC + offset_DAC + 1i*offset_DAC;
end

