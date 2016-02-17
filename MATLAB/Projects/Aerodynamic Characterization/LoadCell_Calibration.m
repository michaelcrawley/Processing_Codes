function [cal,f_indx,T_indx] = LoadCell_Calibration(loadcell)
    switch lower(loadcell)
        case 'jr3_fz_down'
            cal = importdata('JR3_CalibrationMatrix.txt');
            f_indx = @(x) x(:,[3,2,1]).*repmat([-1 1 -1],size(x,1),1);
            T_indx = @(x) x(:,[6,5,4]).*repmat([-1 1 1],size(x,1),1);
        case 'jr3_fz_up'
            cal = importdata('JR3_CalibrationMatrix.txt');
            f_indx = @(x) x(:,[1,2,3]).*repmat([1 -1 1],size(x,1),1);
            T_indx = @(x) x(:,[4,5,6]).*repmat([1 -1 1],size(x,1),1);
        case 'ati_n25_fz_up'
            cal = importdata('ATI_N25_FT14574.txt');
            f_indx = @(x) x(:,[1,2,3]);
            T_indx = @(x) x(:,[4,5,6]);
        case 'ati_n25_fz_down'
            cal = importdata('ATI_N25_FT14574.txt');
            f_indx = @(x) x(:,[1,2,3]).*repmat([1 -1 -1],size(x,1),1);
            T_indx = @(x) x(:,[4,5,6]).*repmat([1 -1 -1],size(x,1),1);
        otherwise
            error('Incorrect load cell definition');
    end
end