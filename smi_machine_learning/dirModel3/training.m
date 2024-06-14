clear
clc

%% 构造数据集
load 'D:\matlab save\self-mixing\smi_machine_learning\dirModel3\dirModel3_small_data_noDC.mat';
X = fringeData';  % 数据集的输入
inputDim = 4;  % 输入维度(没用)
outputDim = 3;  % 输出维度Y（没用）

inputNum = 9656;  % 输入数据总数
inputNum = 4728;  % 输入数据总数
Y = zeros(3,inputNum);  % 数据集输入对应的标签值
for i=1:inputNum
    if(X(31,i)==-1)
        Y(1,i) = 1;
    elseif(X(31,i)==0)
        Y(2,i) = 1;
    else
        Y(3,i) = 1;
    end
end
X = X(1:30,:);
%% 划分数据集
T = randperm(inputNum);  % 将输入数据打散
x_train = X (:,T(1:0.8*inputNum));  % 输入数据的前80%作为训练集
y_train = Y (:,T(1:0.8*inputNum)); 

x_test = X (:,T(0.8*inputNum+1:end));  % 输入数据的前80%作为训练集
y_test = Y (:,T(0.8*inputNum+1:end)); 

%% 数据归一化，注意，数据归一化的要求是每一列代表一个样本！！
% [x_train,xs_train] = mapminmax(X_train,0,1);  % 把训练输入数据归一化到0-1之间，x_train为归一化后的数据，xs_train方便下次同来做同样的归一化处理
% x_test = mapminmax('apply',X_test,xs_train);  % 把测试数据使用 与xs_train同样的归一化方式进行归一化
% 
% [y_train,ys_train] = mapminmax(Y_train,0,1);  % 把输出数据归一化到0-1之间，y_train为归一化后的数据，ys_train主要方便之后反归一化使用（方便看出原数据的变化）

%% 构建神经网络
net = feedforwardnet([10,10,10]);  
net.trainParam.epochs = 1000;  % 训练次数为1000
% net.trainParam.goal = 1e-5; % 训练目标误差小于1e-5时停止训练
net.trainParam.max_fail = 30;  % 训练终止条件之一，默认值为6。即默认情况下会将样本随即分为3类：训练样本，确认样本和测试样本，当确认样本的误差曲线连续6次迭代不在下降。这时训练终止

net.trainParam.lr = 0.1;  % 训练学习速率为0.1
net.trainparam.show=1;  % 显示中间结果的周期

%% 开始训练
net = train(net,x_train,y_train);

%% 测试网络
y_sim = sim(net,x_test);  % 对训练好的神经网络，进行测试集样本验证
% Y_sim = mapminmax('reverse',y_sim,ys_train);  % 对仿真结果进行反归一化处理，便于观测

%% 结果分析
[predict,index1] = max(y_sim);  % 测试数据分类（每行为一类，Max为1即为该类）
[real,index2] = max(y_test);  % 真实数据分类（每行为一类，Max为1即为该类）

X0 = index1==index2;  % 就看预测为1的准确率
X1 = sum(X0~=0);
[c,size] = size(y_sim);

disp('--------------------预测结果-------------------------');
disp(['预测精度为',num2str(X1/size)]);
