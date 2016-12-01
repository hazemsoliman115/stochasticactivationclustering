function rec_power_w = get_path_loss(ap,aa, p_dBm)
    % L  =  128.1 37.6log (R), R in KM
    r = sqrt((ap.x_pos - aa.x_pos)^2 + (ap.y_pos - aa.y_pos)^2);
    p_l = 128.1 + 37.6*log10(r/1000);
    p_l_w = 10^(p_l/10);
    trx_power_w = 10^((p_dBm-30)/10);
    rec_power_w = trx_power_w/p_l_w;