%%INSTRUCTIONS
% This script is executed after load_vtk
% Run this script
% it should create two figures

xpos = VTKData.points(:,1);
ypos = VTKData.points(:,2);
zpos = VTKData.points(:,3);

% for the circle
angles = linspace(0, 2*pi, 500);
radius = 80;
CenterX = 400;
CenterY = 400;
x = radius * cos(angles) + CenterX;
y = radius * sin(angles) + CenterY;


con_all = [];
con_average = [];
con_average = zeros(80,80);
xbin_position = [];
xbin_position = zeros(80,80);
ybin_position = [];
ybin_position = zeros(80,80);
norm_xpos1 = [];
hyp_xpos1 = [];
norm_ypos1 = [];
hyp_ypos1 = [];
for yn = 1:80
    for xn = 1:80
        xbin_position(xn,yn) = 10*(xn-0.5);
        ybin_position(xn,yn) = 10*(yn-0.5);
        for j = 1:numPoints
            if (xn-1)*10<xpos(j) && xpos(j)<=10*xn && (yn-1)*10<ypos(j) && ypos(j)<=10*yn
                con_all(end+1) = VTKData.concentration(j);
            end
        end
        if ~isempty(con_all)
            con_average(xn,yn) = mean(con_all);
        end
        if con_average(xn,yn) < 30.0 && ~isempty(con_all)
            hyp_xpos1(end+1) = xbin_position(xn,yn);
            hyp_ypos1(end+1) = ybin_position(xn,yn);
        end
        if 30 <= con_average(xn,yn) 
            norm_xpos1(end+1) = xbin_position(xn,yn);
            norm_ypos1(end+1) = ybin_position(xn,yn);
        end
        con_all = [];
    end
end

fig1=figure;
plot(hyp_xpos1,hyp_ypos1,'.', Color='Blue', MarkerSize=10);
hold on;
plot(norm_xpos1,norm_ypos1,'.', Color='Red', MarkerSize=10);
axis([0 800 0 800]);
hold on;
set(gca,'FontSize',20,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
Xl=xlabel('x');
Yl=ylabel('y');
plot(x, y, 'k-', 'LineWidth', 2);
saveas(fig1,'average_concentration_1000.png');

pheno_all = [];
pheno_average = [];
pheno_average = zeros(80,80);
xbin_position = [];
xbin_position = zeros(80,80);
ybin_position = [];
xbin_position = zeros(80,80);
norm_xpos = [];
hyp_xpos = [];
norm_ypos = [];
hyp_ypos = [];
nomut_xpos = [];
nomut_ypos = [];
for yn = 1:80
    for xn = 1:80
        xbin_position(xn,yn) = 10*(xn-0.5);
        ybin_position(xn,yn) = 10*(yn-0.5);
        for j = 1:numPoints
            if (xn-1)*10<xpos(j) && xpos(j)<=10*xn && (yn-1)*10<ypos(j) && ypos(j)<=10*yn
                pheno_all(end+1) = VTKData.contpheno(j);
            end
        end
        if ~isempty(pheno_all)
            pheno_average(xn,yn) = mean(pheno_all);
        end
        if pheno_average(xn,yn) < 0.5 && ~isempty(pheno_all)
            hyp_xpos(end+1) = xbin_position(xn,yn);
            hyp_ypos(end+1) = ybin_position(xn,yn);
        end
        if pheno_average(xn,yn) == 0.5 %< 0.55 && 0.45 < pheno_average(xn,yn)
            nomut_xpos(end+1) = xbin_position(xn,yn);
            nomut_ypos(end+1) = ybin_position(xn,yn);
        end
        if 0.5 < pheno_average(xn,yn)
            norm_xpos(end+1) = xbin_position(xn,yn);
            norm_ypos(end+1) = ybin_position(xn,yn);
        end
        pheno_all = [];
    end
end

fig2=figure;
plot(hyp_xpos,hyp_ypos,'.', Color='Blue', MarkerSize=10);
hold on;
plot(norm_xpos,norm_ypos,'.', Color='Red', MarkerSize=10);
plot(nomut_xpos,nomut_ypos,'.', Color='Black', MarkerSize=10);
axis([0 800 0 800]);
hold on;
set(gca,'FontSize',20,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
Xl=xlabel('x');
Yl=ylabel('y');
plot(x, y, 'k-', 'LineWidth', 2);
saveas(fig2,'average_phenotype_1000.png');

