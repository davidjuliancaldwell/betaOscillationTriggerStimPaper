%%
numInt = 7;
chanInt = 14; 

sid = SIDS{numInt};
load(fullfile([sid 'epSTATS-PP-sig-reref-50-new.mat']))



%%
figure
% this sets the figure to be the whole screen
set(gcf, 'Units', 'Inches', 'OuterPosition', [0 0 8 6]);

prettyline(1e3*t,1e6*awins(:, keeps), label(keeps), colors);
xlim(1e3*[min(t) max(t)]);
yl = ylim;
yl(1) = min(-10, max(yl(1),-140*4));
yl(2) = max(10, min(yl(2),100*4));

yl(1) = min(-10, max(yl(1),-340*4));
yl(2) = max(10, min(yl(2),300*4));

ylim(yl);

xlim([-10 60])
ylim([-500 500])

%  vline(1e3*7/efs);
vline(0);
xlabel('time (ms)');
ylabel('ECoG (uV)');
%                 title(sprintf('EP By N_{CT}: %s, %d, {%s}', sid, chan, suffix{typei}))
if (types(typei) == nullType)
    title(sprintf('%s CEPs for Channel %d, Null Condition',sid,chan))
    leg = {'Pre','Post'};
elseif (types(typei) ~= nullType)
    title(sprintf('EP By N_{CT}: %s, %d, {%s}', sid, chan, suffix{typei}))
    title(sprintf('%s CEPs for Channel %d stimuli in {%s}',sid,chan,suffix{typei}))
    leg = {'Pre'};
    for d = 1:length(labelGroupStarts)
        if d == length(labelGroupStarts)
            leg{end+1} = sprintf('%d<=CT', labelGroupStarts(d));
        else
            leg{end+1} = sprintf('%d<=CT<%d', labelGroupStarts(d), labelGroupEnds(d));
        end
    end
end

leg{end+1} = 'Stim Window';
legend(leg, 'location', 'Southeast')
set(gca,'fontsize',18)