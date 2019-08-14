function [] = Gridme( Iin )
%imagesc display of figure with gridlines
%   Detailed explanation goes here

figure;imagesc(Iin);grid on;
ax = gca;
ax.GridColor = 'w';
ax.GridAlpha = 1;
ax.LineWidth = 1;

end

