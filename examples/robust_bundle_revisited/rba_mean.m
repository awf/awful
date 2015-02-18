%% Compute robust mean of data points

clf

%% Create data
pts = [
    randn(16,2) * .02 + au_bsx([.3 .5])
    rand(10,2)*.7+.3
    [randn(10,1) * .02 + .8, rand(10,1)*.07+.3]
    ];
n = size(pts,1);

subplot(121)
hold off
plot(pts(:,1), pts(:,2), '.');
axis([0 1 0 1])
hold on

%%
sty = {'markersize', 7, 'linewidth', 0.25, 'LineSmoothing','on'};

% mean
mu = mean(pts);
plot(mu(1), mu(2), 'ro', sty{:});
% elementwise median
emu = median(pts);
plot(emu(1), emu(2), 'r+', sty{:});

% True L1 minimum.
residuals = @(x) normrows(pts - repmat(x, n, 1));

E_l1 = @(x) sum(residuals(x));
tic
med = fminunc(E_l1, emu);
toc
plot(med(1), med(2), 'rd', sty{:});
drawnow

% Zollhofer minimum
RADIUS = 0.2;
E_psi= @(x) sum(rba_psi(residuals(x)/RADIUS));
tic
rba = fminunc(E_psi, emu);
toc
plot(rba(1), rba(2), 'r*', sty{:});
drawnow

legend('data', 'mean (L2)', 'medmed', 'L1', 'Zollhöfer');

%% 2d plot (L1)
% subplot(122)
% rng = 0:.007:1;
% [xx,yy] = meshgrid(rng);
% zz_l1 = xx;
% for k=1:numel(xx)
%     zz_l1(k) = E_l1([xx(k) yy(k)]);
% end
% zz = zz_l1;
% hold off
% imagesc(rng, rng, -log(zz-min(zz(:))+1e-5))
% axis xy
% colormap hot
% hold on
% plot(mu(1), mu(2), 'wo', sty{:});
% plot(emu(1), emu(2), 'w+', sty{:});
% plot(med(1), med(2), 'wd', sty{:});
% drawnow

% cameratoolbar setmode orbit

%% Robust kernel
subplot(122)
rng = 0:.007:1;
[xx,yy] = meshgrid(rng);
zz_psi = xx;
for k=1:numel(xx)
    zz_psi(k) = E_psi([xx(k) yy(k)]);
end
zz = zz_psi;
hold off
imagesc(rng, rng, -log(zz-min(zz(:))+1e-5))
hold on
contour(rng, rng, zz, 2.^[-4:.2:3 3.1:.1:10], 'c')
colormap hot
axis xy
hold on
plot(pts(:,1), pts(:,2), 'g.');
plot(mu(1), mu(2), 'wo', sty{:});
plot(emu(1), emu(2), 'w+', sty{:});
plot(med(1), med(2), 'wd', sty{:});
plot(rba(1), rba(2), 'w*', sty{:}, 'markersize', 13);
