close all
clear all
format longE

addpath('/path/to/advanpix/toolbox/');

list = [300, 250, 200];
L = length(list);

Marks = {'r.', 'b*', 'kd'};

if length(Marks) < L
  error('Error: Vector of markers is insufficient.');
end % if

iegs = cell(1, L);
legs = cell(1, L);

for idx = 1:L
  n = list(idx);
  nstr = num2str(n);

  fprintf('Computing eigs for n = %3d ... ', n);

  mp.Digits(n);

  m0name = strcat('data/M0_', nstr, '.mat');
  m1name = strcat('data/M1_', nstr, '.mat');
  m2name = strcat('data/M2_', nstr, '.mat');

  M0 = mp.read(m0name);
  M1 = mp.read(m1name);
  M2 = mp.read(m2name);

  [V, e, k] = polyeig(M0, mp('1i')*M1, M2);

  % We save the most accurate result:
  if (n == max(list))
    V0 = V;
    e0 = e;
    k0 = k;
  end % if ()

  iegs{idx} = e;
  legs{idx} = strcat(' ', nstr);

  fprintf('done.\n');

end % for

figure;
for idx = 1:L
    e = iegs{idx};
    plot(real(e), imag(e), Marks{idx}), grid off, hold on
end % for
hold off
l = legend(legs{:});
set(l, 'box', 'off');
set(l, 'location', 'best');
xlabel('$\mathrm{Re}\,\omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
ylabel('$\mathrm{Im}\,\omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
axis([-6 6 -5 1]);
set(gcf, 'color', 'w');