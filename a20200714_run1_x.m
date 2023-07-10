%% analysis for probabilistic teleportation

foldername = '~/teleportation/';
date = '20200714/';
runN = 'run1';
q = 'x';

loadcal = true;
loadsn = true;
loadff = false;
loadtel = true;
loaddn = false;
runanalysis = true;

snruns = 20;
dnruns = 20;
calruns = 20;
telruns = 100;
ffruns = 20;

L = 5e5;

vlog_upper = 7e5;
vlog_lower = 5.8e5;

if loaddn
    
    hetxA = zeros(L*dnruns,1);
    hetpA = zeros(L*dnruns,1);
    homA = zeros(L*dnruns,1);

    for nn = 1:dnruns

        pos = (nn-1)*L+1;

        filename = sprintf('%s%d%s','dark_',nn-1,'_.dat');

        file = fopen([date filename]);
        A = fread(file);
        fclose(file);
        N = length(A);
        
%         skip = 6;
%         ind1 = 1:skip:N;        
    
        ind1 = 1:2:(2*L)-1;
        ind2 = ind1+2*L;
        ind3 = ind2+2*L;
                
        isneg = A(ind1)>127;
        hetx = A(ind1)*256+A(ind1+1)-isneg*256^2;        
        hetxA(pos:pos+L-1) = hetx;

%         isneg = A(ind1+2)>127;
%         hetp = A(ind1+2)*256+A(ind1+3)-isneg*256^2;
        isneg = A(ind2)>127;
        hetp = A(ind2)*256+A(ind2+1)-isneg*256^2;
        hetpA(pos:pos+L-1) = hetp;

%         isneg = A(ind1+4)>127;
%         hom = A(ind1+4)*256+A(ind1+5)-isneg*256^2;
        isneg = A(ind3)>127;  
        hom = A(ind3)*256+A(ind3+1)-isneg*256^2;
        homA(pos:pos+L-1) = hom;
        
        

        fprintf('Loading dark noise: %d/%d.\n',nn,snruns);

    end

    % calculate normalisation
    muD = mean([hetxA hetpA homA]);
    sigD = std([hetxA hetpA homA]);
    
end


%% shot noise

if loadsn

    hetxA = zeros(L*snruns,1);
    hetpA = zeros(L*snruns,1);
    homA = zeros(L*snruns,1);

    for nn = 1:snruns

        pos = (nn-1)*L+1;

        filename = sprintf('%s%s%d%s',runN,'_shot_',nn-1,'_.dat');

        file = fopen([date filename]);
        A = fread(file);
        fclose(file);
        N = length(A);
        
%         skip = 6;
%         ind1 = 1:skip:N;        
    
        ind1 = 1:2:(2*L)-1;
        ind2 = ind1+2*L;
        ind3 = ind2+2*L;
                
        isneg = A(ind1)>127;
        hetx = A(ind1)*256+A(ind1+1)-isneg*256^2;        
        hetxA(pos:pos+L-1) = hetx;

%         isneg = A(ind1+2)>127;
%         hetp = A(ind1+2)*256+A(ind1+3)-isneg*256^2;
        isneg = A(ind2)>127;
        hetp = A(ind2)*256+A(ind2+1)-isneg*256^2;
        hetpA(pos:pos+L-1) = hetp;

%         isneg = A(ind1+4)>127;
%         hom = A(ind1+4)*256+A(ind1+5)-isneg*256^2;
        isneg = A(ind3)>127;  
        hom = A(ind3)*256+A(ind3+1)-isneg*256^2;
        homA(pos:pos+L-1) = hom;
        
        

        fprintf('Loading shot noise: %d/%d.\n',nn,snruns);        
        
    end

    % calculate normalisation
    mu0 = mean([hetxA hetpA homA]);
    sig0 = std([hetxA hetpA homA]);       
    
end

    %% data
if loadcal

        hetxA = zeros(L*calruns,1);
        hetpA = zeros(L*calruns,1);
        homA = zeros(L*calruns,1);
        
        for nn = 1:calruns
            
            pos = (nn-1)*L+1;
            
            filename = sprintf('%s%s%d%s',runN,'_cal_',nn-1,'_.dat');
            
            file = fopen([date filename]);
            A = fread(file);
            fclose(file);
            N = length(A);
            
            %         skip = 6;
            %         ind1 = 1:skip:N;
            
            ind1 = 1:2:(2*L)-1;
            ind2 = ind1+2*L;
            ind3 = ind2+2*L;
            
            isneg = A(ind1)>127;
            hetx = A(ind1)*256+A(ind1+1)-isneg*256^2;
            hetxA(pos:pos+L-1) = hetx;
            
            %         isneg = A(ind1+2)>127;
            %         hetp = A(ind1+2)*256+A(ind1+3)-isneg*256^2;
            isneg = A(ind2)>127;
            hetp = A(ind2)*256+A(ind2+1)-isneg*256^2;
            hetpA(pos:pos+L-1) = hetp;
            
            %         isneg = A(ind1+4)>127;
            %         hom = A(ind1+4)*256+A(ind1+5)-isneg*256^2;
            isneg = A(ind3)>127;
            hom = A(ind3)*256+A(ind3+1)-isneg*256^2;
            homA(pos:pos+L-1) = hom;
            
            fprintf('Loading calibration: %d/%d.\n',nn,calruns);

        end
        
        % normalise
        H = ([hetxA hetpA homA])./sig0;
        
        muC = mean(H);
        sigC = cov(H);
        
        %xin = mux(2)*sqrt(2);
        %vin = sigx(2,2);
                    
end  

if loadff

        hetxA = zeros(L*ffruns,1);
        hetpA = zeros(L*ffruns,1);
        homA = zeros(L*ffruns,1);
        
        for nn = 1:ffruns
            
            pos = (nn-1)*L+1;
            
            filename = sprintf('%s%s%d%s',runN,'_ff_',nn-1,'_.dat');
            
            file = fopen([date filename]);
            A = fread(file);
            fclose(file);
            N = length(A);
            
            %         skip = 6;
            %         ind1 = 1:skip:N;
            
            ind1 = 1:2:(2*L)-1;
            ind2 = ind1+2*L;
            ind3 = ind2+2*L;
            
            isneg = A(ind1)>127;
            hetx = A(ind1)*256+A(ind1+1)-isneg*256^2;
            hetxA(pos:pos+L-1) = hetx;
            
            %         isneg = A(ind1+2)>127;
            %         hetp = A(ind1+2)*256+A(ind1+3)-isneg*256^2;
            isneg = A(ind2)>127;
            hetp = A(ind2)*256+A(ind2+1)-isneg*256^2;
            hetpA(pos:pos+L-1) = hetp;
            
            %         isneg = A(ind1+4)>127;
            %         hom = A(ind1+4)*256+A(ind1+5)-isneg*256^2;
            isneg = A(ind3)>127;
            hom = A(ind3)*256+A(ind3+1)-isneg*256^2;
            homA(pos:pos+L-1) = hom;
            
            fprintf('Loading gain calibration: %d/%d.\n',nn,ffruns);                      
            
            
        end

        H = ([hetxA hetpA homA])./sig0-muC;
        muff = mean(H);
        
        if q=='x'            
            phi = abs(muff(3)/muff(2));
        elseif q=='p'
            phi = abs(muff(3)/muff(1));
        end
              
                    
end 

if loadtel
   fprintf('Loading teleporter data \n')
        hetxA = zeros(L*telruns,1);
        hetpA = zeros(L*telruns,1);
        homA = zeros(L*telruns,1);
        vlog = zeros(telruns,1);
        disc = [];
        
        for nn = 1:telruns
            
            pos = (nn-1)*L+1;
            
            filename = sprintf('%s%s%d%s',runN,'_',nn-1,'_.dat');
            
            file = fopen([date filename]);
            A = fread(file);
            fclose(file);
            N = length(A);
            
            %         skip = 6;
            %         ind1 = 1:skip:N;
            
            ind1 = 1:2:(2*L)-1;
            ind2 = ind1+2*L;
            ind3 = ind2+2*L;
            
            isneg = A(ind1)>127;
            hetx = A(ind1)*256+A(ind1+1)-isneg*256^2;
            hetxA(pos:pos+L-1) = hetx;
            
            %         isneg = A(ind1+2)>127;
            %         hetp = A(ind1+2)*256+A(ind1+3)-isneg*256^2;
            isneg = A(ind2)>127;
            hetp = A(ind2)*256+A(ind2+1)-isneg*256^2;
            hetpA(pos:pos+L-1) = hetp;
            
            %         isneg = A(ind1+4)>127;
            %         hom = A(ind1+4)*256+A(ind1+5)-isneg*256^2;
            isneg = A(ind3)>127;
            hom = A(ind3)*256+A(ind3+1)-isneg*256^2;
            homA(pos:pos+L-1) = hom;

            vlog(nn) = var(hom);

            if and(vlog(nn)<vlog_upper,vlog(nn)>vlog_lower)
                gooddata(pos:pos+L-1) = true(L,1);
            else
                gooddata(pos:pos+L-1) = false(L,1);
            end
            
            %fprintf('Loading teleportation: %d/%d.\n',nn,telruns);                      
            
            
        end
        
        % normalise
        H = ([hetxA hetpA homA])./sig0-muC;

        % % discard bad locking
        H = H(gooddata,:);
        
        mux = mean(H);
        sigx = cov(H);
        
        xout = mux(3);
        vout = sigx(3,3);
        
                    
end 

%%

if runanalysis
    fprintf('Starting postselection \n')
    rep = 1;

    vkeep = zeros(1,rep);
    mukeep = zeros(1,rep);
    fkeep = zeros(1,rep);
    muin = zeros(1,rep);
    muout = zeros(1,rep);

    for mm = 1:rep

        cutoff = 4.7;
        g = 1.54731;
g=1.2;
        len = length(H);

        %K = mvnrnd(zeros(3,1),sigx,1e8);
        %K = H(1:(L*200),:)+muff*mm/20;
       % K = H(1:(L*150),:);
        K=H;

        %succ = exp((0.5*(K(:,1).^2)-cutoff^2)*(1-1/g^2)); 
        succ = exp((0.5*(K(:,1).^2+K(:,2).^2)-cutoff^2)*(1-1/g^2)); % 
        PS = rand(length(K),1) < succ;
        Xwins = sum(PS);
        winrate = Xwins/len;
        Xkeep = K(PS,:);

        % % throw away those above critical point
        % keep = Xkeep(:,3) < 5;
        % Xkeep2 = Xkeep(keep,:);

        sig = cov(Xkeep);
        mu = mean(Xkeep);

        vkeep(mm) = sig(3,3);
        muK = mean(K);
        muin(mm) = -muK(2)*sqrt(2);
        muout(mm) = mu(3);
        fkeep(mm) = fidelity(muin(mm),muout(mm),sig(3,3))

        
        cv(mm)=sig(3,3)-sig(3,2)*sig(2,3)/sig(2,2);
        transfcoeff(mm)=mu(3)*mu(3)/sig(3,3)/(mux(2)*mux(2));
        transfcoeffalt(mm)=mu(3)/sig(3,3)/(sqrt(2)*mux(2));


    end

    % clf
    % figure(1)
    % histogram(Xkeep(:,2),100,'Normalization','pdf')
    % figure(2)
    % % %histogram(Xkeep(:,3),97,'Normalization','pdf')
    % histfit(Xkeep(:,3),94,'Normal')


    disp('------------------------')
    sigx
    sig
    mux
    mu
    winrate
    Xwins
    fidelity(abs(mux(2)*sqrt(2)),abs(mu(3)),sig(3,3))
    disp('------------------------')



 cv
transfcoeff
transfcoeffalt

return

    % xin = -mux(1)*sqrt(2);
    % xout = mux(3);

    % fidelity(1,1,1.2)

end

clf
x = -5:0.01:5;
edges=[-5.25:0.5:5.25];
plot(x,normpdf(x),'k--','Linewidth',1); hold on;
histogram(Xkeep(:,3),edges,'Normalization','pdf','EdgeAlpha',1,'FaceAlpha',0,'Linewidth',1)
[N,edgesused]=histcounts(Xkeep(:,3),edges,'Normalization','pdf')

[muhat, sighat] = normfit(Xkeep(:,3));
plot(x,normpdf(x,muhat,sighat),'r-','Linewidth',2);


xticks([-4:1:4]);
yticks([0:0.1:0.4]);
%xticklabels({})
%yticklabels({})

xlabel('$x$','Interpreter','latex')
ylabel('\langle x \mid \rho \mid x \rangle')

text(2,0.35,'$g = 1.55$','Fontsize',16,'Interpreter','latex')

box on
plot(x,normpdf(x,0,3),'b-','Linewidth',2);
ylim([0,max(normpdf(x))+0.02])
axis([-4.5 4.5 0 0.4])

ylim([0,max(normpdf(x))+0.02])
pbaspect([1.618 1 1])

set(gca,'Fontsize',16)
set(gcf,'Color',[1 1 1])



% sig = zeros(4);
% sig(1,1) = mean([sigx(1,1) sigp(1,1)]);
% sig(2,2) = mean([sigx(2,2) sigp(2,2)]);
% sig(3,3) = sigp(3,3);
% sig(4,4) = sigx(3,3);
% sig(1,3) = sigp(1,3);
% sig(2,4) = sigx(2,3);
% sig(3,1) = sig(1,3);
% sig(4,2) = sig(2,4);
% 
% covmat = infer_4x4(sig)
% 
% %state = beamsplitter([zeros(4,1) covmat],[1 2],0.5)
% %entanglement([zeros(4,1) covmat],[1 2],'EoF/adesso')
% 
% r = 1;
% c = cosh(2*r);
% s = sinh(2*r);
% 
% state = [zeros(4,1) eye(4)];
% state = squeezer(state,1,[r 0]);
% state = squeezer(state,2,[r pi]);
% 
% 
% state = loss(state,1,[0.7 1]);
% state = loss(state,2,[0.7 1]);
% state = beamsplitter(state,[1 2],0.5)
%   
% %entanglement([zeros(4,1) c0],[1 2],'EoF/adesso')
% 
% % %% theory
% % clf
% % 
% % gff = mean(vertcat(gdata,gdata2));
% % xu = -mean(muin);
% % xin = xu*sqrt(2);
% % glist = linspace(0,3.75,100);
% % 
% % 
% % subplot(3,1,1)
% % plot(gff,-muout,'kx')
% % hold on
% % plot(glist,glist*xu,'k-')
% % 
% % axis([0 4 0 3])
% % 
% % 
% % xlabel('Feedforward gain')
% % ylabel('Output mean')
% % 
% % set(gca,'FontSize',16)
% % grid minor
% % pbaspect([1.618 1 1])
% % 
% % 
% % 
% % 
% % subplot(3,1,2)
% % plot(gff,vout,'kx')
% % hold on
% % plot(glist,1+glist.^2,'k-')
% % 
% % axis([0 4 0 15])
% % 
% % xlabel('Feedforward gain')
% % ylabel('Output variance')
% % 
% % set(gca,'FontSize',16)
% % grid minor
% % pbaspect([1.618 1 1])
% % 
% % 
% % 
% % subplot(3,1,3)
% % plot(gff,c,'kx')
% % hold on
% % plot(glist,glist,'k-')
% % 
% % axis([0 4 0 3.5])
% % 
% % xlabel('Feedforward gain')
% % ylabel('HET-HOM covariance')
% % 
% % set(gca,'FontSize',16)
% % grid minor
% % pbaspect([1.618 1 1])
% % 

