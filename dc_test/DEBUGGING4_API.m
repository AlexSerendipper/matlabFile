% %% 该程序用于调试，其实也是非常好用且重要的
% %% 瞎搞，二次包络提取的结果来用，这个是目前最有盼头的
% 
% %% 全局变量
% clc;
% clear all;
% close all;
% fs = 200000;  % 采样率，即fs(s)采一个点。
% N = 4000;  
% fv = 100;  % 震动频
% alpha = 4;
% 
% %% 产生自混合信号
% figure(1);
% subplot(8,1,1);
% t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
% lambda = 650e-9;  % 波长
% A = 2 * lambda;  % 幅值（✔）
% L0 = 20 * lambda;  % 外腔距离（✔） 
% Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
% beta = 1;  % the amplitude of selfmixing signal
% phi0 = 4*pi*(L0+Lt)/lambda;
% 
% p = zeros(1,N);
% C = [3,3];
%  % C的变化是一个正弦曲线，不能随机数！
% C_lower = C(1);
% C_upper = C(2);
% % 这个乘和加保证了c的上下限，这里可以设置变换的周期！！但是这个变换周期需要长一点否则会报错！！
% x = linspace(0, 3*pi, N);
% c = (C_upper-C_lower)/2 * cos(x) + (C_upper - (C_upper-C_lower)/2);
% % plot(x,c);
% 
% for i = 1:N 
%     C = c(i);
%     p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
% end
% % p = awgn(p,30);  % 10db，加高斯噪声
% % p = p .* (1+0.3*cos(2*pi*100*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
% plot(p);
% p_init = p;
% hold on;
% 
% title(['自混合信号,C=',num2str(C), 'alpha=',num2str(alpha)]);
% 
% %% 找到所有的峰谷值（包含跳变点）
% figure(1);
% [top_p,loc_p] = findpeaks(p);
% % loc_p_convert = N2T(loc_p,N,T);
% [top_v,loc_v] = findpeaks(-p);
% top_v = -top_v;
% % loc_v_convert = N2T(loc_v,N,T);
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);
% 
% %% 包络确定方向！
% subplot(8,1,2);
% diffp = diff(p);  % diff是相邻两个数的差分，对第一个位置补0
% diffp = [0,diffp];
% diff_acosp = diff(acos(p));
% diff_acosp = [0,diff_acosp];
% 
% plot(diffp);
% hold on;
% [top_diffp_p,loc_diffp_p] = findpeaks(diffp);  % 拿到极值和索引值
% [top_diffp_v,loc_diffp_v] = findpeaks(-diffp);
% scatter(loc_diffp_p,top_diffp_p);
% scatter(loc_diffp_v,-top_diffp_v);
% en_top = interp1(loc_diffp_p,top_diffp_p,1:N,'spline');  % 三次样条插值，曲线更平滑
% en_bottom = interp1(loc_diffp_v,-top_diffp_v,1:N,'spline');
% 
% en_median = (en_top - (en_top - en_bottom)/2);
% dir = -sign(en_median);  % 让dir暂时指明方向
% 
% plot(en_top);
% plot(en_bottom);
% plot(dir);
% title("diffp及其包络确定方向")
% 
% 
% %% 利用方向信息，求出最合适的找翻转点的范围，大大增加了鲁棒性
% subplot(8,1,3);
% direction_seg1 = [];  % 方向发生变化的点(ˇ∀ˇ)
% for i = 1:length(dir)-1
%     if dir(i) ~=  dir(i+1)
%         direction_seg1 = [direction_seg1, i];
%     end
% end
% 
% direction_seg2 = direction_seg1 + 1;
% direction_seg = [1, sort([direction_seg1,direction_seg2]), length(dir)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)
% top_diffp_seg = [];
% loc_diffp_seg = [];  % 存储翻转点的区间
% 
% for i = 1 : 2 : length(direction_seg)
%     if dir(direction_seg(i):direction_seg(i+1)) < 0  % 在方向小于0的时候求diffp【两端】的谷值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
%         [tem1, tem2] = findpeaks(-diffp(direction_seg(i):direction_seg(i+1)));
%         top_diffp_seg = [top_diffp_seg, -tem1(1), -tem1(end)];
%         loc_diffp_seg = [loc_diffp_seg, tem2(1)+ direction_seg(i) - 1 , tem2(end) + direction_seg(i) - 1];   % 这里返回的索引是个难点！由于我是分段进行峰值寻找，索引值应当加上初始值（假设在[x,y]中找极值，返回索引是2，第2是极值点,实际上对应的索引值为x+2-1）
%         
%     else
%         [tem3, tem4] = findpeaks(diffp(direction_seg(i):direction_seg(i+1)));  % 在方向大于0的时候求diffp【两端】的峰值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
%         top_diffp_seg = [top_diffp_seg, tem3(1), tem3(end)];
%         loc_diffp_seg = [loc_diffp_seg, tem4(1) + direction_seg(i) - 1, tem4(end) + direction_seg(i) - 1];
%     end
% end
% % plot(diffp);
% top_diffp_seg([1,end]) = [];
% loc_diffp_seg([1,end]) = [];  % 最适合找翻转点的区间
% plot(p);
% hold on;
% scatter(direction_seg,0,"r");
% title("根据方向dir，找出用来求diffp极值的区间（红点）（恒定正负1）")
% 
% subplot(8,1,4);
% plot(diffp);
% hold on;
% scatter(loc_diffp_seg,0,"g");
% title("diffp，在direction>0的时候求极小值！！！，direction<0的时候求极大值！！！这样求出的区间(绿点)内不再包含条纹~~~")
% 
% %% 取出翻转点，与上述方法相同！(分成上翻转点和下翻转点 最好是！)
% subplot(8,1,5);
% top_r = [];
% loc_r = [];
% top_r_upper = [];
% loc_r_upper = [];
% top_r_lower = [];
% loc_r_lower = [];
% for i = 1:2:length(loc_diffp_seg)
%     [temp1,temp2] = findpeaks(p(loc_diffp_seg(i):loc_diffp_seg(i+1)));  % 这里求出的就是上翻转点
%     [temp3,temp4] = findpeaks(-p(loc_diffp_seg(i):loc_diffp_seg(i+1)));  % 这里求出的就是下翻转点
%     
%     top_r_upper = [top_r_upper, temp1];
%     loc_r_upper = [loc_r_upper, temp2 + loc_diffp_seg(i) - 1];
%     
%     top_r_lower = [top_r_lower, -temp3];
%     loc_r_lower = [loc_r_lower, temp4 + loc_diffp_seg(i) - 1];
%     
%     top_r = [top_r, temp1, -temp3];
%     loc_r = [loc_r, temp2 + loc_diffp_seg(i) - 1, temp4 + loc_diffp_seg(i) - 1 ];
% end
% 
% plot(p);
% hold on;
% scatter(loc_r,top_r);
% title("在上述区间内求出翻转点")
% 
% %% 挖去峰值中的跳变点，返回，峰值、谷值、跳变点！
% for i = 1:length(loc_p)
%     for j = 1:length(loc_r)
%         if loc_p(i) == loc_r(j)
%             loc_p(i) = nan;
%             top_p(i) = nan;
%         end
%     end
% end
% 
% % 挖去谷值中的跳变点
% for i = 1:length(loc_v)
%     for j = 1:length(loc_r)
%         if loc_v(i) == loc_r(j)
%             loc_v(i) = nan;
%             top_v(i) = nan;
%         end
%     end
% end
% 
% % 删除数组中的nan
% loc_p(isnan(loc_p))=[];
% loc_v(isnan(loc_v))=[];
% top_p(isnan(top_p))=[];
% top_v(isnan(top_v))=[];
% 
% 
% %% 根据求出的翻转点修正一下dir信息！
% subplot(8,1,6);
% direction = zeros(1,N);
% direction(1) = dir(1);
% k = 1;
% for i = 1:N
%     direction(i) = k;
%     if ismember(i, loc_r)
%         k = k * -1;
%     end
% end
% 
% plot(p);
% hold on;
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);
% scatter(loc_r,top_r);
% plot(direction);
% title("得到纯净的峰谷值和完全正确的方向")
% % direction = -direction;  % 方向
% %% 二次包络提取去直流
% figure(2);
% subplot(6, 1, 1);
% plot(p);
% title("smi and envelop")
% hold on;
% scatter(loc_v,top_v);
% en_top = interp1([loc_p,loc_r_upper],[top_p,top_r_upper],1:N,'spline');  % 插值，曲线更平滑
% % en_bottom = interp1([loc_v,loc_r_lower],[top_v,top_r_lower],1:N,'spline');  % 插值，曲线更平滑
% en_bottom = interp1(loc_v,top_v,1:N,'spline');
% en_median = (en_top - (en_top - en_bottom)/2); 
% plot(en_top,"b");
% plot(en_bottom,"g");
% plot(en_median,"r");
% 
% % 第一次包络提取结果，实际上是把翘起来的给往下拉
% subplot(6, 1, 2);
% p1 = p - en_median; 
% plot(p1);
% title('第一次包络提取结果');
% hold on;
% 
% % 第一一次包络提取结果，我再拉一次！！！！！！！！！！！！！！！！！！
% % subplot(6, 1, 3);
% % [top1,loc1] = findpeaks(p);
% % [top2,loc2] = findpeaks(-p);
% % plot(p1);
% % hold on;
% % scatter(loc2,-top2);
% % 
% % en_top = interp1(loc1,top1,1:N,'spline');  % 插值，曲线更平滑
% % en_bottom = interp1(loc2,-top2,1:N,'spline');
% % en_median = (en_top - (en_top - en_bottom)/2); 
% % plot(en_top,"b");
% % plot(en_bottom,"g");
% % title("smi and envelop")
% % subplot(6, 1, 4);
% % p2 = p - en_median; 
% % plot(p2);
% % title('第一一次包络提取结果');
% % hold on;
% 
% % % 第二次包络提取，实际上是往下拉的部分短了一点，通过除法变长
% % subplot(6, 1, 5);
% % [top3,loc3] = findpeaks(p1);
% % [top4,loc4] = findpeaks(-p1);
% % 
% % en_top2 = interp1(loc3,top3,1:N,'spline');
% % en_bottom2 = interp1(loc4,-top4,1:N,'spline');
% % plot(p1);
% % hold on;
% % scatter(loc3,top3);
% % scatter(loc4,-top4);
% % plot(en_top2);
% % plot(en_bottom2);
% % title("第二次包络提取")
% % 
% % % 第二次包络提取结果
% % subplot(6, 1, 6);
% % p2 = p1 ./ (en_top2);
% % % p2 = 2 * (p2 - min(p2))./(max(p2)- min(p2)) - 1;
% % plot(p2);
% % title("第二次包络提取结果");
% % hold on;
% p = p1; 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% 全局变量
% windowLength = 128; % 窗长
% overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
% window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
% fftLength = 5*windowLength;  % 每个时刻傅里叶变换的长度
% V = 0.65; % 抑制因子2
% 
% %% padding,逆变换不可避免的会让信号变短，要padding
% p_beforepadding = p;
% judge = true; 
% padding = windowLength;  % 逆变换不可避免的会让信号变短，要padding
% while judge
%     if mod(((2 * padding + N - overlapLength) / (windowLength - overlapLength)), 1) == 0
%         judge = false;
%     else
%         padding = padding - 1;
%     end
% end
% p = wextend('1D','sym',p,padding);  % 数据延拓
% % p = wextend('1D','zpd',p,padding);  % 0数据延拓
% % p = [zeros(1,padding) p zeros(1,padding)];
% 
% 
% %% 时频分析（抑制）
% figure(3);
% [TF,F,T] = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （TF每一列就是每一个时间维度的频率）
% weight1 = abs(TF)./max(abs(TF));  % 要把这个矩阵当作权值使用，所以按列先归一化
% 
% % weight1(find(weight1 == 1)) = weight1(find(weight1 == 1)) - 1e-11;  % 算法1：是找到weight1中等于1的位置（理论上每列都有），并都减去num极小值，因为这个后边要用做分母，所以减去一个无穷小量，不要让他产生无穷大
% % weight2 = (weight1.^64)./(1 - weight1.^64); % 算法1：构建一个函数，让原来频率分量大的更大，小的更小，windowLenght此处为抑制因子
% 
% weight2 = weight1;
% weight2( find(weight2 < V) ) = 0; % 算法2：嘎嘎好用
% weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化
% TF_curb = TF .* weight2;  % 获得抑制后的信号
% subplot(7,2,[2,4,6]);
% mesh(abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
% % view(0,90); % 设置初始视角为俯视角
% % mesh(angle(TF));
% title('抑制前');
% 
% %% 拉平时频谱去除包络
% % TF_curb(find(TF_curb~=0)) = TF_curb(find(TF_curb~=0))./ abs(TF_curb(find(TF_curb~=0))) ; % 找到其不为0处的索引,把这些处的值设置除其各自的模值，这样这些复数的模长都为1，只不过相角不同
% % TF_curb(find(abs(TF_curb)>0.65)) = TF_curb(find(abs(TF_curb)>0.65))./ abs(TF_curb(find(abs(TF_curb)>0.65))) ; % 找到其不为0处的索引,把这些处的值设置除其各自的模值，这样这些复数的模长都为1，只不过相角不同
% subplot(7,2,[10,12,14]);
% mesh(abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
% % mesh(angle(TF_curb));
% % view(0,90);
% title('抑制后');
% 
% 
% %% 
% p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
% % p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;  % 葛兄归一化
% p_init2 = p;  % 用p_init2存储未归一化的信号
% 
% p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding
% % p_init2 = p_init2(padding+25:end-padding+24);  % ✔去掉之前的padding!!!!!!!!!!!!!!!!!!!!!!!!!!!这里是最大的误差
% % p_init2 = p_init2(padding+45:end-padding+44);  % ✔去掉之前的padding!!!!!!!!!!!!!!!!!!!!!!!!!!!这里是最大的误差
% % p_init2 = p_init2(padding+25:end-padding+24);  % ✔去掉之前的padding!!!!!!!!!!!!!!!!!!!!!!!!!!!这里是最大的误差
% % p_init2 = p_init2(padding-25:end-padding-26);  % ✔去掉之前的padding!!!!!!!!!!!!!!!!!!!!!!!!!!!这里是最大的误差
% % p_init2 = sgolayfilt(p_init2,4,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 
% 
% 
% subplot(7,2,3);
% plot(p_init2);
% title('未归一化');
% 
% %% 希尔伯特变换重构
% 
% 
% 
% subplot(7, 2, 9);
% inverse_hb = (imag(hilbert(p_init2))) .* direction;  % 得到sin ，使用未归一化的信号，精度更高
% phiF_wrapped = atan(inverse_hb./p_init2);  %  使用未归一化的信号，精度更高
% 
% % inverse_hb = (imag(hilbert(p))) .* direction;  % 使用归一化的信号
% % phiF_wrapped = atan(inverse_hb./p);  %  使用归一化的信号
% 
% plot(phiF_wrapped);
% hold on;
% title("phiF-wrapped")
% 
% %% arctan相位解包裹
% for i = 2:N
%     if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
%         phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
%     elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，加pi
%         phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
%     end
% end
% 
% phiF_reconstruct = phiF_wrapped;
% % subplot(7, 1, 6);
% % plot(phiF_reconstruct)
% % title("重构后的phiF")
% 
% %% 计算一下抑制后，C的值
% [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct); 
% 
% %% 重构
% subplot(7, 2, 11);
% phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha));
% Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
% 
% 
% Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% % Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
% plot(Lt,'k');
% hold on;
% 
% plot(Lt_reconstruct,'r'); ylabel('Disp(um)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
% title(['重构后的信号，C_reconstruct=', num2str(C_reconstruct)]);
% 
% %% 误差分析
% subplot(7,2,13);
% plot((Lt-Lt_reconstruct)*10^9); ylabel('error(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
% RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
% title(['绝对误差，RMSE=', num2str(RMSE)]);
% 
% 
% %%  
% figure(5);
% subplot(4,1,1);
% plot(p_init,'r');
% hold on;
% plot(p_init2,'k');
% plot(direction);
% title("抑制且去padding  后 与原信号对比,红前黑后");
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% just for test 
% 
% % figure(4);
% % subplot(4,1,1);
% % plot(p_beforepadding,'r');
% % hold on;
% % plot(p_afterpadding,'k');
% % title("padding前后对比,红前黑后")
% % 
% % subplot(4,1,2);
% % plot(p_beforeinhibit,'r');
% % hold on;
% % plot(p_afterinhibit,'k')
% % title("抑制前后对比(都带padding),红前黑后");
% % 
% % subplot(4,1,3);
% % plot(p_afterinhibit,'r');
% % hold on;
% % plot(p_init2,'k')
% % title("抑制后去padding对比,红前黑后");
% %  
% % subplot(4,1,4);
% % plot(p_init2,'r');
% % hold on;
% % plot(p_beforepadding,'k')
% % plot(direction);
% % title("抑制且去padding后与原信号对比,红前黑后");
% %  
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% subfunction （selfmixing-power）
% function phiF = solve_phiF(C,phi0,alpha)  % 求解出每一个phi0对应的phiF
%     if C<=1  % 每个phi0的解区间
%         [phiF_min, phiF_max] = bounds1(C,phi0);
%     else
%         [phiF_min, phiF_max] = bounds2(C,phi0,alpha);
%     end
%     
%     excessphaze_equation = @(phiF)phiF-phi0+C*sin(phiF+atan(alpha));
%     
%     if (excessphaze_equation (phiF_min)>0)  % 文章的解释为,phiF_min值可能刚好比零点大一点点点，这时候取phiF_min为近似零点
%         phiF = phiF_min;
%     elseif (excessphaze_equation (phiF_max)<0)
%     	phiF = phiF_max;
%     else
%     	phiF = fzero(excessphaze_equation,[phiF_min,phiF_max]);
%     end  
% end
% 
% %---C < 1时的解区间函数------------------------------------------
% function [phiF_min, phiF_max] = bounds1(C,phi0) 
%     phiF_min = phi0-C;
%     phiF_max = phi0+C;
% end
% 
% %---C > 1时的解区间函数------------------------------------------
% function [phiF_min, phiF_max] = bounds2(C,phi0,alpha)
% persistent m; 
% if isempty (m); m = 0; end
% mlower = ceil ((phi0 + atan (alpha) + acos (1/C)- sqrt (C*C- 1))/(2*pi)- 1.5);
% mupper = floor ((phi0 + atan (alpha)- acos (1/C) + sqrt (C*C- 1))/(2*pi)- 0.5);
% if (m < mlower); m = mlower; end
% if (m > mupper); m = mupper; end
% phiF_min = (2*m+1)*pi + acos (1/C)- atan (alpha); 
% phiF_max = (2*m+3)* pi- acos (1/C)- atan (alpha); 
% end
% 
% 
% %% subfunction2(reconstruction_T_N间转换)
% % 该函数实现，将N转换为对应的t，也就是从T=1，N=10这些简单的时候推出来N(i)与t的关系
% % 其中N为需要转换的点，N为采样点数，T为模拟时间
% function  temp = N2T(temp1,N,T)
% temp = zeros(1,length(temp1));    
%     for i = 1:length(temp1) 
%         temp(i) = T*(temp1(i)-1)/N;
%     end
% end
% 
% function temp = T2N(temp1,N,T)
% temp = zeros(1,length(temp1));    
%     for i = 1:length(temp1) 
%         temp(i) = N*temp1(i)/T + 1;
%     end
% end
% 
% 
% %% subfunction （eliminate_dc,输入一段信号p，及其长度N ，输出去除直流分量后的信号p_）
% function g_ = eliminate_dc(p_)
%     p_= fft(p_);
%     p_(1) = 0;
%     g_ = ifft(p_);
% end