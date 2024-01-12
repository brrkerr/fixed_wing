
Alpha = deg2rad(-40:1:40)';
Beta = deg2rad(-20:5:20)';
da = deg2rad(-45:5:45);
de = deg2rad(-45:5:45);
dr = deg2rad(-45:5:45);

beta_pos = 5; % percent change from default params wrt positive beta
beta_neg = -2; % percent change from default params wrt negative beta
beta_coeff = zeros(length(Beta),1);

alpha_pos = 2; % percent change from default params wrt positive alpha
alpha_neg = -2; % percent change from default params wrt negative alpha
alpha_coeff = zeros(length(Alpha),1);

M = 50;
alpha0 = 0.4712;

sigma_alpha = (1+exp(-M.*(Alpha - alpha0))+exp(M.*(Alpha + alpha0))) ...
    ./((1+exp(-M.*(Alpha-alpha0))).*(1+exp(M.*(Alpha+alpha0))));

for i = 1:length(beta_coeff)
    if i < round(length(beta_coeff)/2)
        beta_coeff(i) = (round(length(beta_coeff)/2) - i)*beta_neg;
    elseif i > round(length(beta_coeff)/2)
        beta_coeff(i) = (i - round(length(beta_coeff)/2))*beta_pos;
    end
end

for i = 1:length(alpha_coeff)
    if i < round(length(alpha_coeff)/2)
        alpha_coeff(i) = (round(length(alpha_coeff)/2) - i)*alpha_neg;
    elseif i > round(length(alpha_coeff)/2)
        alpha_coeff(i) = (i - round(length(alpha_coeff)/2))*alpha_pos;
    end
end

%%

CL_alpha_beta = zeros(length(Alpha), length(Beta));
for i = 1:length(Beta)
    CL_alpha_beta(:,i) = (1-sigma_alpha).*(CL_0+(CL_alpha + CL_alpha*beta_coeff(i)/100).*Alpha) ...
    +sigma_alpha.*(2.*sign(Alpha).*sin(Alpha).^2.*cos(Alpha));
end
figure('name', 'CLab')
hold on
for i = 1:length(Beta)
    plot(rad2deg(Alpha), CL_alpha_beta(:,i))
end
xlabel('Alpha [deg]')
ylabel('CLab')
grid minor

CLp = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CLp(i) = (CL_p + CL_p*alpha_coeff(i)/100);
end
figure('name', 'CLp')
hold on
plot(rad2deg(Alpha), CLp)
xlabel('Alpha [deg]')
ylabel('CLp')
grid minor

CLq = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CLq(i) = (CL_q + CL_q*alpha_coeff(i)/100);
end
figure('name', 'CLq')
hold on
plot(rad2deg(Alpha), CLq)
xlabel('Alpha [deg]')
ylabel('CLq')
grid minor

CLr = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CLr(i) = (CL_r + CL_r*alpha_coeff(i)/100);
end
figure('name', 'CLr')
hold on
plot(rad2deg(Alpha), CLr)
xlabel('Alpha [deg]')
ylabel('CLr')
grid minor

CLda = zeros(length(da), length(Alpha));
for i = 1:length(alpha_coeff)
    CLda(:,i) = (CL_da + CL_da*alpha_coeff(i)/100).*da';
end
figure('name', 'CLda')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(da), CLda(:,i))
end
xlabel('da [deg]')
ylabel('CLda')
grid minor

CLde = zeros(length(de), length(Alpha));
for i = 1:length(alpha_coeff)
    CLde(:,i) = (CL_de + CL_de*alpha_coeff(i)/100).*de';
end
figure('name', 'CLde')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(de), CLde(:,i))
end
xlabel('de [deg]')
ylabel('CLde')
grid minor

CLdr = zeros(length(dr), length(Beta));
for i = 1:length(beta_coeff)
    CLdr(:,i) = (CL_dr + CL_dr*beta_coeff(i)/100).*dr';
end
figure('name', 'CLdr')
hold on
for i = 1:length(Beta)
    plot(rad2deg(dr), CLdr(:,i))
end
xlabel('dr [deg]')
ylabel('CLdr')
grid minor


%%

CD_alpha_beta = zeros(length(Alpha), length(Beta));
for i = 1:length(Beta)
    CD_alpha_beta(:,i) = 0.0437 + ((CL_0 + (CL_alpha/5 + CL_alpha/5*beta_coeff(i)/100).*Alpha).^2)...
        /pi/0.9/wing_span^2/S;
    CD_alpha_beta(:,i) = abs(CD_alpha_beta(:,i)); 
end
figure('name', 'CDab')
hold on
for i = 1:length(Beta)
    plot(rad2deg(Alpha), CD_alpha_beta(:,i))
end
plot(rad2deg(Alpha), CD_0 + CD_alpha.*Alpha)
xlabel('Alpha [deg]')
ylabel('CDab')
grid minor

CDp = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CDp(i) = abs((CD_p + CD_p*alpha_coeff(i)/100));
end
figure('name', 'CDp')
hold on
plot(rad2deg(Alpha), CDp)
xlabel('Alpha [deg]')
ylabel('CDp')
grid minor

CDq = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CDq(i) = abs((CD_q + CD_q*alpha_coeff(i)/100));
end
figure('name', 'CDq')
hold on
plot(rad2deg(Alpha), CDq)
xlabel('Alpha [deg]')
ylabel('CDq')
grid minor

CDr = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CDr(i) = abs((CD_r + CD_r*alpha_coeff(i)/100));
end
figure('name', 'CDr')
hold on
plot(rad2deg(Alpha), CDr)
xlabel('Alpha [deg]')
ylabel('CDr')
grid minor

CDda = zeros(length(da), length(Alpha));
for i = 1:length(alpha_coeff)
    CDda(:,i) = abs((CD_da + CD_da*alpha_coeff(i)/100).*da');
end
figure('name', 'CDda')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(da), CDda(:,i))
end
xlabel('da [deg]')
ylabel('CDda')
grid minor

CDde = zeros(length(de), length(Alpha));
for i = 1:length(alpha_coeff)
    CDde(:,i) = abs((CD_de + CD_de*alpha_coeff(i)/100).*de');
end
figure('name', 'CDde')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(de), CDde(:,i))
end
xlabel('de [deg]')
ylabel('CDde')
grid minor

CDdr = zeros(length(dr), length(Beta));
for i = 1:length(beta_coeff)
    CDdr(:,i) = abs((CD_dr + CD_dr*beta_coeff(i)/100).*dr');
end
figure('name', 'CDdr')
hold on
for i = 1:length(Beta)
    plot(rad2deg(dr), CDdr(:,i))
end
xlabel('dr [deg]')
ylabel('CDdr')
grid minor

%%

Cm_alpha_beta = zeros(length(Alpha), length(Beta));
for i = 1:length(Beta)
    Cm_alpha_beta(:,i) = Cm_0 + (Cm_alpha + Cm_alpha*beta_coeff(i)/100).*Alpha;
end
figure('name', 'Cmab')
hold on
for i = 1:length(Beta)
    plot(rad2deg(Alpha), Cm_alpha_beta(:,i))
end
xlabel('Alpha [deg]')
ylabel('Cmab')
grid minor

Cmp = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Cmp(i) = (Cm_p + Cm_p*alpha_coeff(i)/100);
end
figure('name', 'Cmp')
hold on
plot(rad2deg(Alpha), Cmp)
xlabel('Alpha [deg]')
ylabel('Cmp')
grid minor

Cmq = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Cmq(i) = (Cm_q + Cm_q*alpha_coeff(i)/100);
end
figure('name', 'Cmq')
hold on
plot(rad2deg(Alpha), Cmq)
xlabel('Alpha [deg]')
ylabel('Cmq')
grid minor

Cmr = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Cmr(i) = (Cm_r + Cm_r*alpha_coeff(i)/100);
end
figure('name', 'Cmr')
hold on
plot(rad2deg(Alpha), Cmr)
xlabel('Alpha [deg]')
ylabel('Cmr')
grid minor

Cmda = zeros(length(da), length(Alpha));
for i = 1:length(alpha_coeff)
    Cmda(:,i) = (Cm_da + Cm_da*alpha_coeff(i)/100).*da';
end
figure('name', 'Cmda')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(da), Cmda(:,i))
end
xlabel('da [deg]')
ylabel('Cmda')
grid minor

Cmde = zeros(length(de), length(Alpha));
for i = 1:length(alpha_coeff)
    Cmde(:,i) = (Cm_de + Cm_de*alpha_coeff(i)/100).*de';
end
figure('name', 'Cmde')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(de), Cmde(:,i))
end
xlabel('de [deg]')
ylabel('Cmde')
grid minor

Cmdr = zeros(length(dr), length(Beta));
for i = 1:length(beta_coeff)
    Cmdr(:,i) = (Cm_dr + Cm_dr*beta_coeff(i)/100).*dr';
end
figure('name', 'Cmdr')
hold on
for i = 1:length(Beta)
    plot(rad2deg(dr), Cmdr(:,i))
end
xlabel('dr [deg]')
ylabel('Cmdr')
grid minor

%%

CY_alpha_beta = zeros(length(Alpha), length(Beta));
for i = 1:length(Alpha)
    CY_alpha_beta(i,:) = CY_0 + (CY_beta + CY_beta*alpha_coeff(i)/100).*Beta;
end
figure('name', 'CYab')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(Beta), CY_alpha_beta(i,:))
end
xlabel('Beta [deg]')
ylabel('CYab')
grid minor

CYp = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CYp(i) = (CY_p + CY_p*alpha_coeff(i)/100);
end
figure('name', 'CYp')
hold on
plot(rad2deg(Alpha), CYp)
xlabel('Alpha [deg]')
ylabel('CYp')
grid minor

CYq = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CYq(i) = (CY_q + CY_q*alpha_coeff(i)/100);
end
figure('name', 'CYq')
hold on
plot(rad2deg(Alpha), CYq)
xlabel('Alpha [deg]')
ylabel('CYq')
grid minor

CYr = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    CYr(i) = (CY_r + CY_r*alpha_coeff(i)/100);
end
figure('name', 'CYr')
hold on
plot(rad2deg(Alpha), CYr)
xlabel('Alpha [deg]')
ylabel('CYr')
grid minor

CYda = zeros(length(da), length(Alpha));
for i = 1:length(alpha_coeff)
    CYda(:,i) = (CY_da + CY_da*alpha_coeff(i)/100).*da';
end
figure('name', 'CYda')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(da), CYda(:,i))
end
xlabel('da [deg]')
ylabel('CYda')
grid minor

CYde = zeros(length(de), length(Alpha));
for i = 1:length(alpha_coeff)
    CYde(:,i) = (CY_de + CY_de*alpha_coeff(i)/100).*de';
end
figure('name', 'CYde')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(de), CYde(:,i))
end
xlabel('de [deg]')
ylabel('CYde')
grid minor

CYdr = zeros(length(dr), length(Beta));
for i = 1:length(beta_coeff)
    CYdr(:,i) = (CY_dr + CY_dr*beta_coeff(i)/100).*dr';
end
figure('name', 'CYdr')
hold on
for i = 1:length(Beta)
    plot(rad2deg(dr), CYdr(:,i))
end
xlabel('dr [deg]')
ylabel('CYdr')
grid minor

%%

Cl_alpha_beta = zeros(length(Alpha), length(Beta));
for i = 1:length(Alpha)
    Cl_alpha_beta(i,:) = Cl_0 + (Cl_beta + Cl_beta*alpha_coeff(i)/100).*Beta;
end
figure('name', 'Clab')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(Beta), Cl_alpha_beta(i,:))
end
xlabel('Beta [deg]')
ylabel('Clab')
grid minor

Clp = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Clp(i) = (Cl_p + Cl_p*alpha_coeff(i)/100);
end
figure('name', 'Clp')
hold on
plot(rad2deg(Alpha), Clp)
xlabel('Alpha [deg]')
ylabel('Clp')
grid minor

Clq = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Clq(i) = (Cl_q + Cl_q*alpha_coeff(i)/100);
end
figure('name', 'Clq')
hold on
plot(rad2deg(Alpha), Clq)
xlabel('Alpha [deg]')
ylabel('Clq')
grid minor

Clr = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Clr(i) = (Cl_r + Cl_r*alpha_coeff(i)/100);
end
figure('name', 'Clr')
hold on
plot(rad2deg(Alpha), Clr)
xlabel('Alpha [deg]')
ylabel('Clr')
grid minor

Clda = zeros(length(da), length(Alpha));
for i = 1:length(alpha_coeff)
    Clda(:,i) = (Cl_da + Cl_da*alpha_coeff(i)/100).*da';
end
figure('name', 'Clda')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(da), Clda(:,i))
end
xlabel('da [deg]')
ylabel('Clda')
grid minor

Clde = zeros(length(de), length(Alpha));
for i = 1:length(alpha_coeff)
    Clde(:,i) = (Cl_de + Cl_de*alpha_coeff(i)/100).*de';
end
figure('name', 'Clde')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(de), Clde(:,i))
end
xlabel('de [deg]')
ylabel('Clde')
grid minor

Cldr = zeros(length(dr), length(Beta));
for i = 1:length(beta_coeff)
    Cldr(:,i) = (Cl_dr + Cl_dr*beta_coeff(i)/100).*dr';
end
figure('name', 'Cldr')
hold on
for i = 1:length(Beta)
    plot(rad2deg(dr), Cldr(:,i))
end
xlabel('dr [deg]')
ylabel('Cldr')
grid minor

%%

Cn_alpha_beta = zeros(length(Alpha), length(Beta));
for i = 1:length(Alpha)
    Cn_alpha_beta(i,:) = Cn_0 + (Cn_beta + Cn_beta*alpha_coeff(i)/100).*Beta;
end
figure('name', 'Cnab')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(Beta), Cn_alpha_beta(i,:))
end
xlabel('Beta [deg]')
ylabel('Cnab')
grid minor

Cnp = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Cnp(i) = (Cn_p + Cn_p*alpha_coeff(i)/100);
end
figure('name', 'Cnp')
hold on
plot(rad2deg(Alpha), Cnp)
xlabel('Alpha [deg]')
ylabel('Cnp')
grid minor

Cnq = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Cnq(i) = (Cn_q + Cn_q*alpha_coeff(i)/100);
end
figure('name', 'Cnq')
hold on
plot(rad2deg(Alpha), Cnq)
xlabel('Alpha [deg]')
ylabel('Cnq')
grid minor

Cnr = zeros(length(Alpha),1);
for i = 1:length(alpha_coeff)
    Cnr(i) = (Cn_r + Cn_r*alpha_coeff(i)/100);
end
figure('name', 'Cnr')
hold on
plot(rad2deg(Alpha), Cnr)
xlabel('Alpha [deg]')
ylabel('Cnr')
grid minor

Cnda = zeros(length(da), length(Alpha));
for i = 1:length(alpha_coeff)
    Cnda(:,i) = (Cn_da + Cn_da*alpha_coeff(i)/100).*da';
end
figure('name', 'Cnda')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(da), Cnda(:,i))
end
xlabel('da [deg]')
ylabel('Cnda')
grid minor

Cnde = zeros(length(de), length(Alpha));
for i = 1:length(alpha_coeff)
    Cnde(:,i) = (Cn_de + Cn_de*alpha_coeff(i)/100).*de';
end
figure('name', 'Cnde')
hold on
for i = 1:length(Alpha)
    plot(rad2deg(de), Cnde(:,i))
end
xlabel('de [deg]')
ylabel('Cnde')
grid minor

Cndr = zeros(length(dr), length(Beta));
for i = 1:length(beta_coeff)
    Cndr(:,i) = (Cn_dr + Cn_dr*beta_coeff(i)/100).*dr';
end
figure('name', 'Cndr')
hold on
for i = 1:length(Beta)
    plot(rad2deg(dr), Cndr(:,i))
end
xlabel('dr [deg]')
ylabel('Cndr')
grid minor