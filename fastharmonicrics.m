%{
/*
 * Copyright 2025 Saulo Jorge Beltrao de Queiroz

 * If this work is useful to you, please cite:
 ********* S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete
           fourier transform square index coefficients,” IEEE Signal Process. Mag.
           (Tips & Tricks), 2025, accepted for publication.
********** S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete
           fourier transform rectangular index coefficients,”, avalailable online in https://arxiv.org/abs/2504.12551, 2025.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
%}
pkg load signal %for Octave only

% SIC DFT algorithm
function xhat = compactsic(x)
  sqrtN = sqrt(length(x)); 
  xhat = zeros(1, sqrtN); % Initialize output DFT vector
  for n = 0:sqrtN-1
    xhat(n+1) = 0;
    for r = 0:sqrtN-1
      xhat(n+1) = xhat(n+1) + x(n + 1 + r*sqrtN);
    end
  end
end



% RIC DFT algorithm. Assuming length(x) is power of two
function xhat = compactric(x, C)
  L=length(x)/C;
  xhat = zeros(1, C); % Initialize output DFT vector
  for c = 0:C-1
    xhat(c+1) = 0;
    for l = 0:L-1
      xhat(c+1) = xhat(c+1) + x(c + 1 + l*C);
    end
  end
end
% version to improve numerical stability (using the Kahan's addition algorithm) 
function xhat = compactricstable(x, C)
  L = length(x) / C;
  xhat = zeros(1, C); % Initialize output DFT vector
  for c = 0:C-1
    sum = 0;
    c1 = c + 1; % Shift for MATLAB indexing (1-based index)
    c_offset = c1 * L; % The offset for the group
    
    % Kahan summation for better accuracy
    c_error = 0;
    for l = 0:L-1
      temp = x(c1 + l*C) - c_error;  % Subtract the previous error
      sum = sum + temp;  % Add the current element
      c_error = (sum - (sum + temp)) - temp;  % Compute new error
    end
    
    xhat(c1) = sum;  % Store the final sum for this group
  end
end

% Use for the computation of a specific frequency
function X_k = compute_dft_frequency(x, N, k)
  X_k = 0;                
  for n = 0:N-1
    X_k = X_k + x(n+1) * exp(-2j * pi * k * n / N);
  end
end

[x, Fs] = audioread('piano_A4_Fs44100Hz.wav'); %fs 7744
q = 13;
Nfft=2^q;
%=================================================================>
for p = 9:12
    Nsic = 7744;    % number of points to analyze

    C=2^p;
    L=2^(q-p);
    sqrtN=sqrt(Nsic);

    alpha = 0.5; %for tukey windowing
    %windowfft = tukeywin(Nfft, alpha); 
    windowfft = hann(Nfft); 
    x_windowedfft = x(1:Nfft) .* windowfft';
    tic
    yfft=fft(x_windowedfft);
    toc

%    xhatric=compactricstable(x(1:Nfft), C); 
    xhatric=compactric(x(1:Nfft), C); 
%    windowric = tukeywin(C, alpha);
    windowric = hann(C);
    xhatwindowedric = xhatric(1:C) .* windowric';
    Yric=zeros(1,Nfft);
    Yric(Yric == 0) = 10e-8;
    Yricaux=fft(xhatwindowedric);

    for c = 1:C
        Yric((c-1)*L+1) = Yricaux(c);
    end

    xhat=compactsic(x(1:Nsic)); 
    windowsic = tukeywin(sqrtN, alpha);
    xhatwindowed = xhat(1:sqrtN) .* windowsic';

    Y=zeros(1,Nsic);
    Y(Y == 0) = 10e-5;

    Y(1)=compute_dft_frequency(xhatwindowed, sqrtN, 0); %DC frequency
    nHarmonics=6;
    for k = 1:nHarmonics
        Y(k*sqrtN + 1) = compute_dft_frequency(xhatwindowed, sqrtN, k); 
    end

    deltaFfft = Fs / Nfft;
    deltaFric = Fs / Nfft;
    deltaFsic = Fs / Nsic;

    frequenciessic = (0:Nsic/2-1) * deltaFsic;
    frequenciesric = (0:Nfft/2-1) * deltaFric;
    frequenciesfft = (0:Nfft/2-1) * deltaFfft;
    magnitudefft = abs(yfft(1:Nfft/2)) / Nfft;
    magnituderic = abs(Yric(1:Nfft/2)) / Nfft;
    magnitudesic = abs(Y(1:Nsic/2)) / Nsic;

    figure;
    semilogy(frequenciesfft(1:Nfft/2), magnitudefft(1:Nfft/2),'LineWidth', 1.2);  
    hold on
    semilogy(frequenciesric(1:Nfft/2), magnituderic(1:Nfft/2));  
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');

    %=====> Dynamic legend with C value
    ric_label = sprintf('%d-point RIC', C);
    l = legend('8192-point FFT', ric_label);
    set(l, 'FontSize', 14);
    set(gca, 'FontSize', 14);

    grid on;
    axis([0 4400 10^-5 10])
    xticks(0:440:4400)

    filename = sprintf('A4pianokey%d.eps', C);
    print(filename, '-depsc');
end

