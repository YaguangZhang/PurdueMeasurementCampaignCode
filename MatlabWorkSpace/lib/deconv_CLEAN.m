function h_t = deconv_CLEAN(x_t, y_t, cutoff_lev_dB, iter_num_max)

% h_t = deconv_td_subtr(x_t, y_t, cutoff_lev_dB, iter_num_max)
%
% This function performs subtractive deconvolution
%
% h_t is the estimate of the impulse response such that y_t = conv(x_t, h_t)
% -freq_cutoff_lev_dB is the level below which the correlation is considered negligible

x_t_len = length(x_t);

auto_corr_x = conv(fliplr(x_t), x_t);
norm_factor = max(auto_corr_x);

cross_corr_x_y = conv(fliplr(x_t), y_t);

samp_lev = 10^(-cutoff_lev_dB/20)*max(abs(cross_corr_x_y));

h_t = zeros(1, length(y_t) - x_t_len + 1);

contin_deconv = 1;
iter_num = 0;
while (contin_deconv)
    iter_num = iter_num + 1;
    
    peak_index = min(find(abs(cross_corr_x_y(x_t_len:end-x_t_len+1)) == max(abs(cross_corr_x_y(x_t_len:end-x_t_len+1)))));
    peak = cross_corr_x_y(peak_index + x_t_len - 1);
    
    h_t(peak_index) = h_t(peak_index) + 1/norm_factor*peak;
    
    cross_corr_x_y(peak_index:peak_index+2*x_t_len-2) = cross_corr_x_y(peak_index:peak_index+2*x_t_len-2) - peak/norm_factor*auto_corr_x;
   
    if ((max(abs(cross_corr_x_y(x_t_len:end-x_t_len+1))) < samp_lev) | (iter_num == iter_num_max))
        contin_deconv = 0;
    end;    
end;

% iter_num

if (iter_num == iter_num_max)
    disp('Error: The algorithm could not complete in the maximum number of iterations specified.');
    h_t = [];
end;