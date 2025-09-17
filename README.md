# week3-Z-transform

clear; clc; close all;
saveFigs = true;         % set false to skip saving PNGs
figFmt   = '-dpng';      % PNG format for exportgraphics/print

%% A1 — Finite Sequences → Polynomials (and verify)
disp('=== A1: Finite sequences as polynomials in z^{-1} ===');
syms n z

% i) x[n] = {1,2,5} at n={0,1,2}
X1 = 1 + 2*z^(-1) + 5*z^(-2);   % explicit polynomial in z^{-1}
disp('X1(z) ='); pretty(X1)

% ii) x[n] = {0,3,0,4} at n={0,1,2,3}
X2 = 0 + 3*z^(-1) + 0*z^(-2) + 4*z^(-3);
disp('X2(z) ='); pretty(X2)

% ROC note for finite right-sided sequences (printed to console)
disp('A1 ROC note: Finite-length right-sided sequences have ROC = entire z-plane except possibly z=0.');

% (Optional illustration: build from deltas)
% x[n] = sum_k x[k] z^{-k}  (already reflected in X1, X2)

%% A2 — Infinite Sequences & ROC
disp('=== A2: Infinite sequences & ROC ===');

% (a) x[n] = a^n u[n], a = 0.6  -> X(z) = 1/(1 - a z^{-1}), ROC |z| > |a|
a = sym(0.6);
Xa = ztrans(a^n * heaviside(n), n, z);
Xa_s = simplify(Xa);
disp('(a) a=0.6, X(z)='); pretty(Xa_s);
disp('ROC: |z| > 0.6');

% (b) x[n] = (-0.8)^n u[n] -> X(z) = 1/(1 + 0.8 z^{-1}), ROC |z| > 0.8
b = sym(-0.8);
Xb = ztrans(b^n * heaviside(n), n, z);
Xb_s = simplify(Xb);
disp('(b) a=-0.8, X(z)='); pretty(Xb_s);
disp('ROC: |z| > 0.8');

% (c) x[n] = -(0.9)^n u[-n-1] (left-sided) -> X(z) = z/(z-0.9), ROC |z| < 0.9
c = sym(0.9);
Xc = ztrans( - (c^n) * heaviside(-n-1), n, z);
Xc_s = simplify(Xc);
disp('(c) left-sided, X(z)='); pretty(Xc_s);
disp('ROC: |z| < 0.9');

disp('ROC & poles: Right-sided -> ROC outside outermost pole; Left-sided -> ROC inside innermost pole; ROC excludes poles.');

%% A3 — Linearity & Shifting
disp('=== A3: Linearity & Time-Shift properties ===');

x1 = (sym(0.5))^n * heaviside(n);
x2 = (sym(-0.5))^n * heaviside(n);

X1 = simplify(ztrans(x1,n,z));    % 1/(1 - 0.5 z^{-1})
X2 = simplify(ztrans(x2,n,z));    % 1/(1 + 0.5 z^{-1})
X_lin = simplify(ztrans( 2*x1 - 3*x2, n, z));
X_shift = simplify(ztrans( subs(x1,n,n-3), n, z));  % Z{x1[n-3]} = z^{-3} X1(z)

disp('X1(z)='); pretty(X1);
disp('X2(z)='); pretty(X2);
disp('Z{2x1[n] - 3x2[n]}='); pretty(X_lin);
disp('Z{x1[n-3]}='); pretty(X_shift);
disp('ROC (for both): |z| > 0.5 (right-sided).');

%% A4 — Inverse Z-Transform
disp('=== A4: Inverse Z-transform ===');

Xa = 1 / (1 - 0.7*z^(-1));
xa = simplify(iztrans(Xa, z, n));
disp('(a) xa[n] ='); pretty(xa);
disp('By inspection: pole at 0.7, ROC |z|>0.7 -> xa[n]=0.7^n u[n].');

Xb = (1 - 0.5*z^(-1)) / (1 - 0.8*z^(-1));
xb = simplify(iztrans(Xb, z, n));
disp('(b) xb[n] ='); pretty(xb);
disp('By long division/partial fraction: xb[n]=0.625*delta[n] + 0.375*(0.8)^n u[n].');

%% A5 — H(z), Poles/Zeros & Frequency Response
disp('=== A5: H(z), poles/zeros, frequency response ===');

b = [1 -2.4  2.88];               % numerator coefficients
a = [1 -0.8  0.64];               % denominator coefficients

% (a) Numeric poles & zeros
zeros_H = roots(b);
poles_H = roots(a);
disp('Zeros (numeric):'); disp(zeros_H.');
disp('Poles (numeric):'); disp(poles_H.');

% Pole-zero plot
f1 = figure('Name','A5_pz'); zplane(b,a); grid on; title('Pole-Zero Plot of H(z)');
if saveFigs, exportgraphics(f1, 'A5_pz.png'); end

% (b) Magnitude & phase via freqz
[H,w] = freqz(b,a,512,'whole');             % whole 0..2pi for completeness
w = w(1:257); H = H(1:257);                 % keep 0..pi for plotting
f2 = figure('Name','A5_freq');
subplot(2,1,1); plot(w/pi, abs(H), 'LineWidth',1); grid on;
ylabel('|H(e^{j\omega})|'); title('Magnitude Response'); xlim([0 1]);
subplot(2,1,2); plot(w/pi, angle(H), 'LineWidth',1); grid on;
xlabel('\omega/\pi'); ylabel('Phase (rad)'); xlim([0 1]);
if saveFigs, exportgraphics(f2, 'A5_freq.png'); end

% (c) Quick time-domain test
n0 = 0:511;
x = sin(0.2*pi*n0) + 0.5*sin(0.8*pi*n0);    % two-tone input
y = filter(b,a,x);

f3 = figure('Name','A5_inout_x');
plot(n0, x, 'LineWidth',1); grid on; xlabel('n'); ylabel('x[n]');
title('Input x[n] = sin(0.2\pi n) + 0.5 sin(0.8\pi n)');
if saveFigs, exportgraphics(f3, 'A5_input.png'); end

f4 = figure('Name','A5_inout_y');
plot(n0, y, 'LineWidth',1); grid on; xlabel('n'); ylabel('y[n]');
title('Output y[n] = filter(b,a,x)');
if saveFigs, exportgraphics(f4, 'A5_output.png'); end

% Short textual classification (also print to console)
disp('Filter type comment:');
disp(['Poles near radius ~0.8 at ~±60° produce a resonant peak around ~0.35π; ', ...
      'response is narrowband/peaking (resonator-like), stable (|p|<1), and non-minimum-phase due to zeros outside unit circle.']);

%% Console summary for README copy-paste
disp('=== Summary (copy for README) ===');
fprintf('A1: X1(z) = 1 + 2 z^{-1} + 5 z^{-2};  X2(z) = 3 z^{-1} + 4 z^{-3}. ROC: entire plane except possibly z=0.\n');
fprintf('A2: a^n u[n] -> 1/(1 - a z^{-1}), ROC |z|>|a|.  (-0.8)^n u[n] -> 1/(1+0.8 z^{-1}), ROC |z|>0.8.\n');
fprintf('    Left-sided -(0.9)^n u[-n-1] -> z/(z-0.9), ROC |z|<0.9.\n');
fprintf('A3: Z{2x1 - 3x2} and Z{x1[n-3]} verified symbolically; ROC |z|>0.5.\n');
fprintf('A4: 1/(1-0.7 z^{-1}) -> 0.7^n u[n].  (1-0.5 z^{-1})/(1-0.8 z^{-1}) -> 0.625 δ[n] + 0.375 (0.8)^n u[n].\n');
fprintf('A5: Poles: '); disp(poles_H.'); fprintf('Zeros: '); disp(zeros_H.');
fprintf('Saved figures: A5_pz.png, A5_freq.png, A5_input.png, A5_output.png\n');

![image](https://github.com/Khan548-codes/week3-Z-transform/blob/main/images/s1.png)
![image](https://github.com/Khan548-codes/week3-Z-transform/blob/main/images/s2.png)
![image](https://github.com/Khan548-codes/week3-Z-transform/blob/main/images/s3.png)
![image](https://github.com/Khan548-codes/week3-Z-transform/blob/main/images/s4.png)

# A1 – Finite Sequences to Polynomials

Plots: There are no real plots here, just polynomials. Each sequence turns into a polynomial in 
𝑧^−1 where the numbers in the sequence are just the coefficients.

Reflection: I learned that any finite sequence can be written as a short polynomial in 
𝑧^−1
. It was interesting to see how each term in the sequence matches directly with a power of z^−1
. The ROC is simple for these because the sequences end after a few samples.



# A2 – Infinite Sequences and ROC

Plots: If we plotted them, the first one would look like a decaying exponential, the second one like a decaying exponential that flips sign every step, and the third one would be a left-sided sequence that grows backwards.

Reflection: I noticed the ROC is very important. For right-sided signals the ROC is outside the pole, but for left-sided ones it is inside. It made sense after checking the definitions. Also, the alternating sign in the second case makes the plot go up and down.




# A3 – Linearity and Shifting
•	# Plots: If plotted, the exponentials just combine into another signal, and shifting makes the signal start later (it is zero until the shift).
•	Reflection: I understood how linearity makes things easy: you can just take the Z-transform of each part and add them. Time-shifting is also clear because it just adds a factor of z^-k. I also saw that the ROC stays the same when signals are shifted.





# A4 – Inverse Z-Transform
•	Plots: The first inverse gives a decaying exponential, while the second one has a spike at n=0 and then an exponential tail.
•	Reflection: I learned how to go back from X(z)to x[n]. The method with partial fractions or inspection shows where the impulse or exponential terms come from. It was nice to see how the poles decide the behavior of the sequence.




# A5 – H(z), Poles/Zeros and Frequency Response
•	Plots:
o	The pole-zero plot shows poles inside the unit circle and zeros outside.
o	The magnitude plot has a peak around 0.35π.
o	The phase plot changes smoothly with frequency.
o	The input signal has two sinusoids, but the output mainly keeps the low-frequency part because of the filter.
•	Reflection: I noticed how the poles close to the unit circle make a resonance in the frequency response. The filter is stable because the poles are inside the unit circle, but it is not minimum phase because some zeros are outside. The time-domain test showed how the filter really amplifies one frequency and reduces the other.


