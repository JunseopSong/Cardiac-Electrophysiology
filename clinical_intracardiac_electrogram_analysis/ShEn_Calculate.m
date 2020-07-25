clear all;

egm = load('egm2.csv');

for i=1:3000
    u_bin(i) = -15 + 0.01*i;
    y_bin(i) = -15 + 0.01*(i-1);
end

nde = 256;
egmLength = 2035;
ShEn = zeros(nde ,1);

for nodeId = 1:nde
    disp(nodeId);
    
    EGM = squeeze(egm(:,nodeId));
    
    for i=1:3000
        bin(i)=0;
    end
    
    % Search the 1st bin
    i=1;
    while(u_bin(i)<=EGM(1))
        i=i+1;
    end
    ci=i;
    
    xo=0;
    ind=0;
    for i=1:(egmLength-1)  % 1000 ms
        
        slope = EGM(i+1) - EGM(i);
        if( slope>0 )
            
            x = (u_bin(ci)-EGM(i))/slope + (i-1);
            
            while (x<i)
                bin(ci+ind) = bin(ci+ind) + ( x - xo );
                
                xo=x;
                ind = ind+1;
                
                x = (u_bin(ci+ind)-EGM(i))/slope + (i-1);
            end
            
            bin(ci+ind) = bin(ci+ind) + ( i - xo );
            
            xo=i;
            ci=ci+ind;
            if(x==i)
                ci=ci+1;
            end
            ind=0;
        elseif(slope<0)
            x = (y_bin(ci)-EGM(i))/slope + (i-1);
            
            while (x<i)
                bin(ci-ind) = bin(ci-ind) + ( x - xo );
                
                xo=x;
                ind = ind+1;
                
                x = (y_bin(ci-ind)-EGM(i))/slope + (i-1);
            end
            
            bin(ci-ind) = bin(ci-ind) + ( i - xo );
            
            xo=i;
            ci=ci-ind;
            if(x==i)
                ci=ci-1;
            end
            ind=0;
        else
            bin(ci)=bin(ci)+1;
            xo=i;
        end
        
    end
    
    
    sum = egmLength-1;
    ShEn(nodeId,1)=0;
    for i=1:3000
        prob(i)=bin(i)/sum;
        if(prob(i)>0)
            ShEn(nodeId,1) = ShEn(nodeId,1) - prob(i)*log2(prob(i));
        end
    end
    
end


fid = fopen('ShEn2.txt', 'w');
fprintf(fid, '%f\n', ShEn);
fclose(fid);
