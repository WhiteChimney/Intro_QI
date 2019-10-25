% 此文件可显示HWP与QWP作用在量子态上后在Bloch球上的轨迹
% 其中白点为初态，黑点为x,y,z轴正向部分与Bloch球的交点
% 其余点为末态随波片角度变化而形成的轨迹
% 注：三维图可以拖动查看不同视角

%% 预设参数
C1 = 1; C2 = 0; phi_0 = pi*0/4;   % 末态表示为 C1|0> + C2*exp(i*phi_0)|1>
alpha = 0:pi/1000:pi*16/16;         % 波片角度，起始：步距（精度）：结束
HWPorQWP = 1;                       % 选择波片，0为HWP，1为QWP

%% 计算末态及对应Bloch球上的坐标
% 定义HWP与QWP函数及其逆
HWP =@(theta) [cos(2*theta),  sin(2*theta);...
               sin(2*theta), -cos(2*theta)];
QWP =@(theta) [cos(2*theta)-1i,   -sin(2*theta);...
                  sin(2*theta), cos(2*theta)+1i];
HWP_r =@(theta) HWP(theta);
QWP_r =@(theta) 0.5*QWP(theta)';

% 入射态
psi_0 = [C1; C2*exp(1i*phi_0)];
psi_0 = RealUp(psi_0);

% 求经波片后的态在Bloch球上所对应的坐标
xyz_0 = BlochCoord(psi_0);             % 初始态坐标
xyz_x = zeros(length(alpha),3);        % 末态坐标
for i = 1:length(alpha)
    if HWPorQWP == 0
        psi_x = HWP(alpha(i)) * psi_0;
    else
        psi_x = QWP(alpha(i)) * psi_0;
    end
    psi_x = RealUp(psi_x);
    xyz_x(i,:) = BlochCoord(psi_x);
end

% 画图
hold on
[sx,sy,sz] = sphere;                             % Bloch球
surf(sx,sy,sz,'FaceAlpha',0.5)
shading interp
text(xyz_0(1),xyz_0(2),xyz_0(3),'o','color','w') % 白点为初态
text(1,0,0,'o','color','k')                      % 黑点为x,y,z坐标轴正向
text(0,1,0,'o','color','k')
text(0,0,1,'o','color','k')
scatter3(xyz_x(:,1),xyz_x(:,2),xyz_x(:,3),'.')   % 画出末态随波片角度变化的轨迹
hold off


%% 子函数

function res = RealUp(psi)
% 将输入的态正规化，即归一且第一个系数为实数
    res = [psi(1)*conj(psi(1)); psi(2)*conj(psi(1))];
    res = res/norm(res);
end

function xyz = BlochCoord(psi)
% 输入一个态，返回其在Bloch球上的直角坐标
    % Pauli矩阵
    pauli_x = [0,1;1,0];
    pauli_y = [0,-1i;1i,0];
    pauli_z = [1,0;0,-1];
    % 坐标对应公式
    xyz = [psi'*pauli_x*psi;...
           psi'*pauli_y*psi;...
           psi'*pauli_z*psi];
end
