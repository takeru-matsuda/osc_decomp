function ret = plot_phase_area(t,phase,lphase,uphase,C)
    for i=1:length(t)-1
        lturn = 0;
        uturn = 0;
        if (phase(i)-lphase(i))*(phase(i+1)-lphase(i+1))<0
            if lphase(i)<phase(i) && phase(i)-phase(i+1)<pi
                lturn = -1;
            end
            if lphase(i)>phase(i) && phase(i)-phase(i+1)>-pi
                lturn = 1;
            end
        else
            if phase(i)-phase(i+1)>pi
                lturn = 1;
            end
            if phase(i)-phase(i+1)<-pi
                lturn = -1;
            end
        end
        if (phase(i)-uphase(i))*(phase(i+1)-uphase(i+1))<0
            if uphase(i)<phase(i) && phase(i)-phase(i+1)<pi
                uturn = -1;
            end
            if uphase(i)>phase(i) && phase(i)-phase(i+1)>-pi
                uturn = 1;
            end
        else
            if phase(i)-phase(i+1)>pi
                uturn = 1;
            end
            if phase(i)-phase(i+1)<-pi
                uturn = -1;
            end
        end
        if lphase(i)<uphase(i)
            if lphase(i+1)<uphase(i+1)
                ld = lphase(i);
                lu = uphase(i);
                rd = lphase(i+1)+lturn*2*pi;
                ru = uphase(i+1)+uturn*2*pi;
                fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                ld = lphase(i)-lturn*2*pi;
                lu = uphase(i)-uturn*2*pi;
                rd = lphase(i+1);
                ru = uphase(i+1);
                fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
            else
                ld = lphase(i);
                lu = uphase(i);
                rd = lphase(i+1)+lturn*2*pi;
                ru = uphase(i+1)+uturn*2*pi;
                fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                ld = lphase(i)-(lturn+uturn)*2*pi;
                lu = uphase(i)-(lturn+uturn)*2*pi;
                rd = lphase(i+1)-uturn*2*pi;
                ru = uphase(i+1)-lturn*2*pi;
                fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
            end
        else
            if lphase(i+1)<uphase(i+1)
                ld = lphase(i)-lturn*2*pi;
                lu = uphase(i)-uturn*2*pi;
                rd = lphase(i+1);
                ru = uphase(i+1);
                fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                ld = lphase(i)+uturn*2*pi;
                lu = uphase(i)+lturn*2*pi;
                rd = lphase(i+1)+(lturn+uturn)*2*pi;
                ru = uphase(i+1)+(lturn+uturn)*2*pi;
                fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
            else
                if lturn == 0
                    ld = lphase(i);
                    lu = uphase(i)+2*pi;
                    rd = lphase(i+1);
                    ru = uphase(i+1)+2*pi;
                    fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                    ld = lphase(i)-2*pi;
                    lu = uphase(i);
                    rd = lphase(i+1)-2*pi;
                    ru = uphase(i+1);
                    fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                else
                    ld = lphase(i);
                    lu = uphase(i)+lturn*2*pi;
                    rd = lphase(i+1)+lturn*2*pi;
                    ru = uphase(i+1)+lturn*4*pi;
                    fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                    ld = lphase(i)-lturn*2*pi;
                    lu = uphase(i);
                    rd = lphase(i+1);
                    ru = uphase(i+1)+lturn*2*pi;
                    fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                    ld = lphase(i)-lturn*4*pi;
                    lu = uphase(i)-lturn*2*pi;
                    rd = lphase(i+1)-lturn*2*pi;
                    ru = uphase(i+1);
                    fill([t(i) t(i+1) t(i+1) t(i)]',[ld rd ru lu]',C,'EdgeColor','none');
                end
            end
        end
    end
    ylim([-pi pi]);
    ret = 0;
end

