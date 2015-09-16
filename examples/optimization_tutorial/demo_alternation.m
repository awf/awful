function log_data = demo_alternation(x0, y0, xrange, yrange, lz, pmin, maxiter)

if nargin == 0
  %% Test  
  evalin('base', 'demo_alternation(x0, y0, xrange, yrange, z, pmin, 19);');
  return
end


hvlinestyles = {'color', [1 1 1]*0, 'linestyle', '-', 'linewidth', 4};
xk = x0;
yk = y0;
[xx,yy] = meshgrid(xrange, yrange);
log_data = [];
hline = @(y) plot([min(xrange) max(xrange)], [y y]);
vline = @(x) plot([x x], [min(yrange) max(yrange)]);

for iter = 1:maxiter
  clf
  colormap jet
  imagesc(xrange, yrange, lz(xx,yy)); axis xy
  hold on
  err = lz(xk,yk)/lz(pmin(1),pmin(2)); err = (err - 1)*100;
  th = text(5, 5, sprintf('Count = %d\nError = %.1f%%', iter-1, err), 'fontsize', 16);
  plot(xk,yk,'wo', 'linewidth', 4);
  
  if (iter < 10) % || (iter < 40 && (rem(iter, 10) == 0))
    rest = @() pause;
  else
    rest = @() drawnow;
  end

  if rem(iter, 2)
    set(hline(yk), hvlinestyles{:});
    xkp1 = fminbnd(@(x) lz(x, yk), min(xrange), max(xrange));
    rest();
    plot(xkp1, yk, 'wo', 'linewidth', 4)
    xk = xkp1;
  else
    set(vline(xk), hvlinestyles{:});
    ykp1 = fminbnd(@(y) lz(xk, y), min(yrange), max(yrange));
    rest();
    plot(xk, ykp1, 'wo', 'linewidth', 4);
    yk = ykp1;
  end
  log_data = [log_data; lz(xk,yk)];
  err = lz(xk,yk)/lz(pmin(1),pmin(2)); err = (err - 1)*100;
  delete(th);
  text(5, 5, sprintf('Count = %d\nError = %.1f%%', iter, err), 'fontsize', 16);
  %rest();
end
