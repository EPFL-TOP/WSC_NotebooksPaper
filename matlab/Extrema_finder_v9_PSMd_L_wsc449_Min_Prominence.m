% Calculating T_plus and T_minus data


load('/Volumes/orange2/wsc_raw_intensities_matlab_1906/ws_v8_WSC449_RI_BG.mat')


her1 = horzcat(wsc449_RI_PSMd_L_yfp(1:121,:));


%%%%%%% Minimum time distance %%%%%%%%%%%%

Min_peak = 40;


%mesp_onset = [];

[m,n] = size(her1);


time = (1:m-1)'*10;

% The rising time (from minima to maxima), the falling time (opposite),
% ratios for TP and TM, differences (or amplitudes) for Tp and Tm.


T_plus_vec = [];   % This is the vector of Tp
T_minus_vec = [];  % This is the vector of Tm
T_period_vec = []; % This is the vector of Tp + Tm
Time_peak_vec =[];
max_val_vec_1 = [];
max_val_vec_2 = [];
min_val_vec = [];
Int_max_vec = []; % COPY
Int_min_vec = []; % COPY
loc_min_vec = []; %This is the time of the minima
Ratio_plus_vec = []; %This is the ratio corresponding to Tp
Ratio_minus_vec = [];%This is the ratio corresponding to Tm
Diff_plus_vec = [];  %This the differences corresponding to Tp
Diff_minus_vec = []; %This the differences corresponding to Tm
Last_Peak_time_vec = [];  % This stores the last peak time
Last_Peak_value_vec = [];  % This stores the last peak time
Max_Peak_time_vec = []; 
Max_peak_value_vec = [];
Integral_max_peak_vec = [];
Integral_vec = [];
Second_Last_Peak_time_vec = [];  % This stores the last peak time
Second_Last_Peak_value_vec = [];  % This stores the last peak time
Freq_index_vec =[];
Last_Diff_vec = [];  % This stores the last difference
Second_Last_Diff_vec = []; % This stores the second last difference



% Enter the initial time frame where for each her1 signal (time of first
% maxima minus one)

% Starting_point_PSMd_L_yfp = [];
% 
% End_point_PSMd_L_yfp = [];

% load('/Volumes/orange/wsc_raw_intensities_matlab/ws_WSC416_RI_BG.mat')

Starting_point_PSMd_L = vertcat(Starting_point_PSMd_L,Starting_point_PSMd_L);

End_point_PSMd_L = vertcat(End_point_PSMd_L,End_point_PSMd_L);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting_point_PSMd_L_yfp = Starting_point_PSMd_L_yfp-1;
% 
% End_point_PSMd_L_yfp = End_point_PSMd_L_yfp-1;


% This is where all the magic happens (Save a copy of the whole file before modifying)

Cases2 = [];
Cases =[];
l=0;

cmap = jet(n);


for i = 1 : n
    
if sum(ismember(Cases2,i)) == 1
    
else
    
% Interpolation of measurement


[max_val,loc_max] = findpeaks(her1(Starting_point_PSMd_L(i,1):End_point_PSMd_L(i,1),i),'MinPeakProminence',Min_peak);    

[min_val,loc_min] = findpeaks(-her1(Starting_point_PSMd_L(i,1):End_point_PSMd_L(i,1),i),'MinPeakProminence',Min_peak);


%Finding Maxima or minima accurately

[m,n2] = size(loc_max);


for j = 1 : m
    
   c(j,:) = polyfit(time(Starting_point_PSMd_L(i,1)+loc_max(j,1)-2:Starting_point_PSMd_L(i,1)+loc_max(j,1),1),her1(Starting_point_PSMd_L(i,1)+loc_max(j,1)-2:Starting_point_PSMd_L(i,1)+loc_max(j,1),i),2);
    
   t_max(j,1) = -c(j,2)/(2*c(j,1));
   
   Int_max(j,1) = c(j,3) - c(j,2)^2/(4*c(j,1));
      
end

clear c

[m,n2] = size(loc_min);


for j = 1 : m
    
   k(j,:) = polyfit(time(Starting_point_PSMd_L(i,1)+loc_min(j,1)-2:Starting_point_PSMd_L(i,1)+loc_min(j,1),1),-her1(Starting_point_PSMd_L(i,1)+loc_min(j,1)-2:Starting_point_PSMd_L(i,1)+loc_min(j,1),i),2);
    
   t_min(j,1) = -k(j,2)/(2*k(j,1));
   
   Int_min(j,1) = -(k(j,3) - k(j,2)^2/(4*k(j,1)));
   
end


clear k

min_val = -min_val;


Last_Peak_time_vec(i,1) = time(loc_max(end,1))+time(Starting_point_PSMd_L(i,1),1);  
Last_Peak_value_vec(i,1) = max_val(end,1);
Second_Last_Peak_time_vec(i,1) = time(loc_max(end-1,1))+time(Starting_point_PSMd_L(i,1),1);
Second_Last_Peak_value_vec(i,1) = max_val(end-1,1);
Last_Diff_vec(i,1) = max_val(end,1)-min_val(end,1);  % This stores the last difference



[M,I] = max(max_val);
Max_peak_value_vec(i,1) = M;
Max_Peak_time_vec(i,1) = time(loc_max(I,1),1)+time(Starting_point_PSMd_L(i,1),1); 
Integral_max_peak_vec(i,1) = sum(her1(Starting_point_PSMd_L(i,1):loc_max(I,1),i));
Integral_vec(i,1)= sum(her1(Starting_point_PSMd_L(i,1):End_point_PSMd_L(i,1),i));


max_val_vec_1 = vertcat(max_val_vec_1,max_val(1:end-1,1));
max_val_vec_2 = vertcat(max_val_vec_2,max_val(2:end,1));

%Getting a new way to find maxima

Int_max_vec = vertcat(Int_max_vec,Int_max); %%%%% COPY
Int_min_vec = vertcat(Int_min_vec,Int_min);


if sum(ismember(Cases,i)) == 1
    
%     l = l + 1;

min_val_vec = vertcat(min_val_vec,min_val(1:end-1,1));

loc_min_vec = vertcat(loc_min_vec,t_min(:,1));   %CHANGED

Time_peak_vec =vertcat(Time_peak_vec,t_max);

T_minus = t_min(1:end-1,1) - t_max(1:end-1,1);

T_plus = t_max(2:end,1) - t_min(1:end-1,1);

Ratio_plus = Int_max(2:end,1)./Int_min(1:end-1,1);

Ratio_plus2 = Int_max(1:end-1,1)./Int_min(1:end-1,1);

Diff_plus = Int_max(2:end,1)-Int_min(1:end-1,1)+1;

Diff_plus2 = Int_max(1:end-1,1)-Int_min(1:end-1,1)+1;

else
    
    l = l + 1;

   min_val_vec = vertcat(min_val_vec,min_val);
  
   loc_min_vec = vertcat(loc_min_vec,t_min(:,1)); %CHANGED
   
   %Time_peak_vec =vertcat(Time_peak_vec,t_max(2:end,1))
   
   T_minus = t_min(1:end,1) - t_max(1:end-1,1);
   
   T_plus = t_max(2:end,1) - t_min(1:end,1);
   
   Ratio_plus = Int_max(2:end,1)./Int_min(1:end,1);
   
   Ratio_plus2 = Int_max(1:end-1,1)./Int_min(1:end,1);
   
   Diff_plus = Int_max(2:end,1)-Int_min(1:end,1)+1; % Diff_plus = Int_max(1:end-1,1)-Int_min(1:end,1);
   
   Diff_plus2 = Int_max(1:end-1,1)-Int_min(1:end,1)+1;  % Diff_plus2 = Int_max(1:end-1,1)-Int_min(1:end,1);

    
end



T_period = T_minus + T_plus;

Index_period = T_period(2:end)./T_period(1:end-1);

Freq_index_vec =vertcat(Freq_index_vec,Index_period);

T_period_vec = vertcat(T_period_vec,T_period);

T_minus_vec = vertcat(T_minus_vec,T_minus);

T_plus_vec = vertcat(T_plus_vec,T_plus);

Ratio_plus_vec = vertcat(Ratio_plus_vec,Ratio_plus);

Ratio_minus_vec = vertcat(Ratio_minus_vec,Ratio_plus2);

Diff_plus_vec = vertcat(Diff_plus_vec,Diff_plus);
   
Diff_minus_vec = vertcat(Diff_minus_vec,Diff_plus2);

Time_peak_vec =vertcat(Time_peak_vec,t_max);

clear t_max t_min Int_max Int_min Diff_plus Diff_plus2

[Size_T(l,1),r] = size(T_minus);

end



end






% These are the plots you get (you can get more!)


figure
subplot(1,2,1)
plot(Diff_plus_vec,T_plus_vec,'xk','MarkerSize',8,'LineWidth',1)
xlabel('x^+-x^-','fontsize',18)
ylabel('T^+','fontsize',18)
set(gca,'FontSize',18)
ylim([0 100])
ylim([0 100])
subplot(1,2,2)
plot(Diff_minus_vec,T_minus_vec,'xk','MarkerSize',8,'LineWidth',1)
xlabel('x^+-x^-','fontsize',18)
ylabel('T^-','fontsize',18)
set(gca,'FontSize',18)
ylim([0 100])

%Generating the PDFs

% [bandwidth1,density1,X1,Y1]=kde2d([Diff_plus_vec/100,T_plus_vec/10]);
% [bandwidth2,density2,X2,Y2]=kde2d([Diff_minus_vec/100,T_minus_vec/10]);
% 
% 
% figure
% subplot(1,2,1)
% imagesc(X1(1,65:174),Y1(:,1),density1(:,65:174))
% set(gca, 'YDir', 'normal')
% subplot(1,2,2)
% imagesc(X2(1,65:174),Y2(:,1),density2(:,65:174))
% set(gca, 'YDir', 'normal')
% 
% 
% Rho1 = sum(density1);
% Rho1 = repmat(Rho1,256,1);
% Rho2 = sum(density2);
% Rho2 = repmat(Rho2,256,1);
% 
% Cond_up = density1(:,65:174)./Rho1(:,65:174);
% Cond_down = density2(:,65:174)./Rho2(:,65:174);
% 
% Norm_up = repmat(sum(Cond_up)*(Y1(2,1)-Y1(1,1)),256,1);
% Norm_down = repmat(sum(Cond_down)*(Y2(2,1)-Y2(1,1)),256,1);
% 
% Cond_up = Cond_up./Norm_up;
% Cond_down = Cond_down./Norm_down;
% 
% mean_up = sum(Y1(:,1).*Cond_up)*(Y1(2,1)-Y1(1,1));
% mean_down = sum(Y2(:,1).*Cond_down)*(Y2(2,1)-Y2(1,1));
% 
% var_up = sum(Y1(:,1).^2.*Cond_up)*(Y1(2,1)-Y1(1,1))-mean_up.^2;
% var_down = sum(Y2(:,1).^2.*Cond_down)*(Y2(2,1)-Y2(1,1))-mean_down.^2;
% 
% 
% figure
% subplot(2,2,1)
% imagesc(X1(1,65:174),Y1(:,1),Cond_up)
% set(gca, 'YDir', 'normal')
% subplot(2,2,2)
% imagesc(X2(1,65:174),Y2(:,1),Cond_down)
% set(gca, 'YDir', 'normal')
% subplot(2,2,3)
% plot(X1(1,65:174),mean_up,X2(1,65:174),mean_down)
% subplot(2,2,4)
% plot(X2(1,65:174),var_up,X2(1,65:174),var_down)
% 
