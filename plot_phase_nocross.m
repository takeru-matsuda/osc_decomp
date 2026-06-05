function plot_phase_nocross(t,phase,lineoption,linewidth)
    N = length(t);
    
    t_extended = zeros(1, 5*N);
    phase_extended = zeros(1, 5*N);
    
    idx = 1;
    t_extended(idx) = t(1);
    phase_extended(idx) = phase(1);
    
    for i=1:N-1
        if phase(i)-phase(i+1) > pi
            turn = t(i) + (pi-phase(i))/(phase(i+1)+2*pi-phase(i))*(t(i+1)-t(i));
            t_extended(idx+1:idx+4) = [turn, NaN, turn, t(i+1)];
            phase_extended(idx+1:idx+4) = [pi, NaN, -pi, phase(i+1)];
            idx = idx + 4;
        elseif phase(i)-phase(i+1) < -pi
            turn = t(i) + (-pi-phase(i))/(phase(i+1)-2*pi-phase(i))*(t(i+1)-t(i));
            t_extended(idx+1:idx+4) = [turn, NaN, turn, t(i+1)];
            phase_extended(idx+1:idx+4) = [-pi, NaN, pi, phase(i+1)];
            idx = idx + 4;
        else
            t_extended(idx+1) = t(i+1);
            phase_extended(idx+1) = phase(i+1);
            idx = idx + 1;
        end
    end
    
    t_plot = t_extended(1:idx);
    phase_plot = phase_extended(1:idx);
    
    plot(t_plot, phase_plot, lineoption, 'LineWidth', linewidth);
end
