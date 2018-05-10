%              _       _   _
%   __  _____ | | ___ | |_| |
%   \ \/ / _ \| |/ _ \| __| |
%    >  < (_) | | (_) | |_| |
%   /_/\_\___/|_|\___/ \__|_|
%
% help: shows the activation functions of channel
%
function ax = show(conductance,ax)


[m_inf, h_inf, tau_m, tau_h] = xolotl.getGatingFunctions(conductance);

V = linspace(-100,100,1e3);

% evaluate these functions 
minf = NaN*V;
hinf = NaN*V;
taum = NaN*V;
tauh = NaN*V;

Ca = 3e3;

for i = 1:length(V)
	if nargin(m_inf) == 1
		minf(i) = m_inf(V(i));
	else
		minf(i) = m_inf(V(i),Ca);
	end
	if nargin(h_inf) == 1
		hinf(i) = h_inf(V(i));
	else
		hinf(i) = h_inf(V(i),Ca);
	end
	
	taum(i) = tau_m(V(i));
	tauh(i) = tau_h(V(i));
end

if nargin < 3


	% check if there is any figure that we can plot on
	allfigs = get(0, 'Children'); 
	f = [];
	for i = 1:length(allfigs)
		if ~isempty(allfigs(i).Tag) && strcmp(allfigs(i).Tag,'xolotl-show')
			f = allfigs(i);
			ax  = findall(f,'type','axes');

		end
	end

	if isempty(f)

		f = figure('outerposition',[100 100 1000 900],'PaperUnits','points','PaperSize',[1000 500]); hold on
		for i = 1:4
			ax(i) = subplot(2,2,i); hold on
		end
		ax(1).Tag = 'm_inf';
		ax(2).Tag = 'h_inf';
		ax(3).Tag = 'tau_m';
		ax(4).Tag = 'tau_h';

		ylabel(ax(1),'m_{inf}')
		xlabel(ax(1),'V (mV)')

		xlabel(ax(2),'V (mV)')
		ylabel(ax(2),'h_{inf}')

		ylabel(ax(3),'\tau_{m} (ms)')
		xlabel(ax(3),'V (mV)')
		set(ax(3),'YScale','log')

		ylabel(ax(4),'\tau_{h} (ms)')
		xlabel(ax(4),'V (mV)')
		set(ax(4),'YScale','log')

		f.Tag = 'xolotl-show';
	end
end


plot(ax(find(strcmp({ax.Tag},'m_inf'))),V,minf,'DisplayName',conductance);
plot(ax(find(strcmp({ax.Tag},'h_inf'))),V,hinf,'DisplayName',conductance);
plot(ax(find(strcmp({ax.Tag},'tau_m'))),V,taum,'DisplayName',conductance);
plot(ax(find(strcmp({ax.Tag},'tau_h'))),V,tauh,'DisplayName',conductance);


prettyFig();

axes(ax(1))
legend;

% turn all YLim modes to auto
for i = 1:length(ax)
	ax(i).YLimMode = 'auto';
end

