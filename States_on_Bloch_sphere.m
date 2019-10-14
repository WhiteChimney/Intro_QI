% ���ļ�����ʾHWP��QWP����������̬�Ϻ���Bloch���ϵĹ켣
% ���а׵�Ϊ��̬���ڵ�Ϊx,y,z�����򲿷���Bloch��Ľ���
% �����Ϊĩ̬�沨Ƭ�Ƕȱ仯���γɵĹ켣
% ע����άͼ�����϶��鿴��ͬ�ӽ�

%% Ԥ�����
C1 = 1; C2 = 0; phi_0 = pi*0/4;   % ĩ̬��ʾΪ C1|0> + C2*exp(i*phi_0)|1>
alpha = 0:pi/1000:pi*16/16;         % ��Ƭ�Ƕȣ���ʼ�����ࣨ���ȣ�������
HWPorQWP = 1;                       % ѡ��Ƭ��0ΪHWP��1ΪQWP

%% ����ĩ̬����ӦBloch���ϵ�����
% ����HWP��QWP����������
HWP =@(theta) [cos(2*theta),  sin(2*theta);...
               sin(2*theta), -cos(2*theta)];
QWP =@(theta) [cos(2*theta)-1i,   -sin(2*theta);...
                  sin(2*theta), cos(2*theta)+1i];
HWP_r =@(theta) [cos(2*theta),  sin(2*theta);...
                 sin(2*theta), -cos(2*theta)];
QWP_r =@(theta) 0.5*[cos(2*theta)+1i,    sin(2*theta);...
                       -sin(2*theta), cos(2*theta)-1i];

% ����̬
psi_0 = [C1; C2*exp(1i*phi_0)];
psi_0 = RealUp(psi_0);

% �󾭲�Ƭ���̬��Bloch��������Ӧ������
xyz_0 = BlochCoord(psi_0);             % ��ʼ̬����
xyz_x = zeros(length(alpha),3);        % ĩ̬����
for i = 1:length(alpha)
    if HWPorQWP == 0
        psi_x = HWP(alpha(i)) * psi_0;
    else
        psi_x = QWP(alpha(i)) * psi_0;
    end
    psi_x = RealUp(psi_x);
    xyz_x(i,:) = BlochCoord(psi_x);
end

% ��ͼ
hold on
[sx,sy,sz] = sphere;                             % Bloch��
surf(sx,sy,sz,'FaceAlpha',0.5)
shading interp
text(xyz_0(1),xyz_0(2),xyz_0(3),'o','color','w') % �׵�Ϊ��̬
text(1,0,0,'o','color','k')                      % �ڵ�Ϊx,y,z����������
text(0,1,0,'o','color','k')
text(0,0,1,'o','color','k')
scatter3(xyz_x(:,1),xyz_x(:,2),xyz_x(:,3),'.')   % ����ĩ̬�沨Ƭ�Ƕȱ仯�Ĺ켣
hold off


%% �Ӻ���

function res = RealUp(psi)
% �������̬���滯������һ�ҵ�һ��ϵ��Ϊʵ��
    res = [psi(1)*conj(psi(1)); psi(2)*conj(psi(1))];
    res = res/norm(res);
end

function xyz = BlochCoord(psi)
% ����һ��̬����������Bloch���ϵ�ֱ������
    rela_amp = psi(2)/psi(1);
    theta = 2*atan(abs(rela_amp));
    phi = mod(imag(log(rela_amp/abs(rela_amp))),2*pi);
    xyz = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
end
