%% 用于手动确认方向信息，输入翻转点arr数组(顺序输入) ,返回自混合信号方向（方向有误乘-1即可）
function direction = SMI_API_DIR(arr,N)
    arr1 = [1, sort([arr, arr+1]), N];
    temp = zeros(1,N);
    flag = 1;
    
    for i = 1:2:length(arr1)-1      
        temp(arr1(i):arr1(i+1)) = flag;
        flag = flag * -1;
    end
    direction = sign(temp);
end

