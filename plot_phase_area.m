function plot_phase_area(t,lphase,uphase,color)

    N = length(t);
    hold on;

    for i = 1:N-1
        ti  = t(i); 
        ti1 = t(i+1);
        
        lp0 = lphase(i);  lp1 = lphase(i+1);
        up0 = uphase(i);  up1 = uphase(i+1);
        
        while up0 < lp0, up0 = up0 + 2*pi; end
        while up0 - lp0 > 2*pi, up0 = up0 - 2*pi; end
        
        while lp1 < lp0 - pi,   lp1 = lp1 + 2*pi; end
        while lp1 > lp0 + pi,   lp1 = lp1 - 2*pi; end
        
        while up1 < up0 - pi,   up1 = up1 + 2*pi; end
        while up1 > up0 + pi,   up1 = up1 - 2*pi; end
        
        while up1 < lp1, up1 = up1 + 2*pi; end
        while up1 - lp1 > 2*pi, up1 = up1 - 2*pi; end
        
        x_poly = [ti; ti1; ti1; ti];
        y_poly = [lp0; lp1; up1; up0];
        
        for shift = -2:2
            y_shift = y_poly + shift * 2 * pi;
            
            if any(y_shift >= -pi & y_shift <= pi) || ...
               (min(y_shift) < -pi && max(y_shift) > pi)
                
                y_clipped = max(-pi, min(pi, y_shift));
                
                if max(y_clipped) - min(y_clipped) > 1e-5
                    fill(x_poly, y_clipped, color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
                end
            end
        end
    end
    
    ylim([-pi pi]);
end
