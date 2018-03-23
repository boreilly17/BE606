function [CancerElementConcentrations] = p53model(t, x, c)

p53g = x(1);
p53m = x(2);
p53p2 = x(3);
p53p4 = x(4);
mdm2ga = x(5);
mdm2gb = x(6);
mdm2ma = x(7);
mdm2p = x(8);
mdm2mb = x(9);
mdm2p_mdm2p = x(10);
mdm2p_p53p2 = x(11);
mdm2p_mdm2p_p53p2 = x(12);
p53p2u = x(13);

d_p53g = 0;
d_p53m = (c(1) * p53g) - (c(5) * p53m);
d_p53p2 = (c(2) * p53m) + (2*c(4) * p53p4) - (c(3) * p53p2 * p53p2) - ...
    (c(6) * p53p2) - (c(29) * mdm2p * p53p2) + (c(30) * mdm2p_p53p2) + ...
    (c(47) * p53p2) - (c(35) * mdm2p_mdm2p * p53p2) - ...
    (c(36) * mdm2p_mdm2p_p53p2);
d_p53p4 = (c(3) * p53p2 * p53p2) - (c(7) * p53p4) - (c(4) * p53p4);
d_mdm2ga = (c(11) * mdm2gb) - (c(10) * mdm2ga);
d_mdm2gb = (c(10) * mdm2ga) - (c(11) * mdm2gb);
% Note: following equation has delay1 on first term
d_mdm2ma = (c(8) * mdm2ga) - (c(14) * mdm2ma);
d_mdm2p = (c(9) * mdm2ma) - (c(16) * mdm2p) + 2*(c(18) * mdm2p_mdm2p) - ...
    (c(17) * mdm2p * mdm2p) - (c(29) * mdm2p * p53p2) + ...
    (c(30) *mdm2p_p53p2) + (c(13) * mdm2mb);
%Note: following equation has delay2 on first term
d_mdm2mb = (c(12) * mdm2gb) - (c(15) * mdm2mb);
d_mdm2p_mdm2p = (c(17) * mdm2p * mdm2p) - (c(18) * mdm2p_mdm2p)...
    - (c(19) * mdm2p_mdm2p) - ...
    (c(35) * mdm2p_mdm2p * p53p2) + (c(36) * mdm2p_mdm2p_p53p2) + ...
    (c(45) * mdm2p_mdm2p_p53p2);
d_mdm2p_p53p2 = (c(29) * mdm2p * p53p2) - (c(30) * mdm2p_p53p2) - ...
    (c(31) * mdm2p_p53p2 * mdm2p) + (c(32) * mdm2p_mdm2p_p53p2);
d_mdm2p_mdm2p_p53p2 = (c(31) * mdm2p_p53p2 * mdm2p) - ...
    (c(32) * mdm2p_mdm2p_p53p2) + (c(35) * mdm2p_mdm2p * p53p2) - ...
    (c(36) * mdm2p_mdm2p_p53p2) - (c(45) * mdm2p_mdm2p_p53p2);
d_p53p2u = (c(45) * mdm2p_mdm2p_p53p2) - (c(47) * p53p2u) - ...
    (c(48) * p53p2u);

CancerElementConcentrations = [d_p53g;
                                d_p53m;
                                d_p53p2;
                                d_p53p4;
                                d_mdm2ga;
                                d_mdm2gb;
                                d_mdm2ma;
                                d_mdm2p;
                                d_mdm2mb;
                                d_mdm2p_mdm2p;
                                d_mdm2p_p53p2;
                                d_mdm2p_mdm2p_p53p2;
                                d_p53p2u;
                                ];

end
