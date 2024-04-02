%% 输入大概的电压，实现自动定标
%% 注意输入的fv和fs需要准确，其他的如phi和offset的自动定标范围是可以手动修改的
function [Lt_standard,mv_standard,phi_standard,offset_standard] = SMI_API_STANDARD_AUTO(Lt_reconstruct,mv, fv, fs, N) 
    RMSE = 99;
    Lt_standard = [];
    mv_standard = [];
    phi_standard = [];
    offset_standard = [];
    for i = mv-10:mv+10
        for phi=1:128
            for offset = -20:0.1:20
                 [Lt] = MOVE_API_STANDARD(mv, fv, fs, N, phi, offset);
                 temp_rmse = sqrt(mean((Lt-Lt_reconstruct).^2));
                 if temp_rmse<RMSE
                     RMSE = temp_rmse;
                     Lt_standard = Lt;
                     mv_standard = mv;
                     phi_standard = phi;
                     offset_standard = offset;
                 end
            end
        end
    end
end

