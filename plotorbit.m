function [lineobj] = plotorbit(ax, r0, v0, dt, res, name, color)
tpl = linspace(0, dt, res);
rpl = zeros(res, 3);
for index=1:res
    rpl(index,:) = kepler_univar(r0, v0, tpl(index), 1)';
    if index > 1
        if any(abs(rpl(index,:)./rpl(index-1,:)-1)>0.1)
            rpl(index,:) = kepler_univar(r0, v0, tpl(index), 1, false, 0.5)';
        end
    end
end
lineobj = plot3(ax, rpl(:,1), rpl(:,2), rpl(:,3), 'LineStyle','-','MarkerSize',10,'Color',color,'DisplayName',name);
end
