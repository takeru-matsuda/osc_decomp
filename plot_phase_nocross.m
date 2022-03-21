function ret = plot_phase_nocross(t,phase,lineoption,linewidth)
    for i=1:length(t)-1
        if phase(i)-phase(i+1) > pi
            turn = t(i)+(pi-phase(i))/(phase(i+1)+2*pi-phase(i))*(t(i+1)-t(i));
            plot([t(i) turn],[phase(i) pi],lineoption,'LineWidth',linewidth);
            plot([turn t(i+1)],[-pi phase(i+1)],lineoption,'LineWidth',linewidth);
        else
            if phase(i)-phase(i+1) < -pi
                turn = t(i)+(phase(i)+pi)/(phase(i)-phase(i+1)+2*pi)*(t(i+1)-t(i));
                plot([t(i) turn],[phase(i) -pi],lineoption,'LineWidth',linewidth);
                plot([turn t(i+1)],[pi phase(i+1)],lineoption,'LineWidth',linewidth);
            else
                plot([t(i) t(i+1)],[phase(i) phase(i+1)],lineoption,'LineWidth',linewidth);
            end
        end
    end
    ret = 0;
end

