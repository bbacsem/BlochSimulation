% 영상 크기 (세로: 101, 가로: 11)
ny = 101;
nx = 101;

% [-1,1] 범위 좌표 생성 (phantom()와 유사한 좌표계)
x_lin = linspace(-1, 1, nx);
y_lin = linspace(-1, 1, ny);
[X, Y] = meshgrid(x_lin, y_lin);

%% Tissue 영역 마스크 생성
% Skull: 전체 머리 외곽 (a=b=1)에서 내부의 뇌 영역(a=b=0.85) 제외
skull_mask = ellipse_mask(X, Y, 0, 0, 1.0, 1.0, 0) & ~ellipse_mask(X, Y, 0, 0, 0.85, 0.85, 0);

% Grey Matter: skull 내부 영역(a=b=0.85)에서 중앙 white matter(a=b=0.6) 제외
grey_mask  = ellipse_mask(X, Y, 0, 0, 0.85, 0.85, 0) & ~ellipse_mask(X, Y, 0, 0, 0.6, 0.6, 0);

% White Matter: 중앙에 위치한 영역 (a=b=0.6)
white_mask = ellipse_mask(X, Y, 0, 0, 0.6, 0.6, 0);

% Ventricle: white matter 내부의 뇌실 (예시: 오른쪽 약간 치우친 위치)
vent_mask = ellipse_mask(X, Y, 0.3, 0.1, 0.2, 0.3, 0);
% white matter 영역에서 뇌실 영역은 제외
white_mask = white_mask & ~vent_mask;

%% 채널별 intensity 설정
% T1 영상: 일반적으로 CSF(ventricle)는 낮은 신호, white matter가 가장 밝음
T1 = 0.3 * skull_mask + 1.6 * grey_mask + 1.3 * white_mask + 4 * vent_mask;

% T2 영상: CSF(ventricle)는 높은 신호, white matter는 낮은 신호
T2 = 0.8e-3 * skull_mask + 45e-3 * grey_mask + 40e-3 * white_mask + 90e-3 * vent_mask;

% Proton Density (PD): white와 grey matter는 비교적 밝게, CSF도 밝은 편
PD = 0.2 * skull_mask + 0.8 * grey_mask + 0.9 * white_mask + 0.95 * vent_mask;

% 3채널 이미지 생성 (101×11×3)
phantom_img = cat(3, T1, T2, PD);
phantom_img(phantom_img==0) = NaN;

phantom_T1 = phantom_img(:,:,1);
phantom_T2 = phantom_img(:,:,2);
phantom_PD = phantom_img(:,:,3);
save('phantom101.mat','phantom_T1','phantom_T2','phantom_PD');
%% 결과 display
figure;
subplot(2,3,1);
imshow(T1, []);
title('T1');

subplot(2,3,2);
imshow(T2, []);
title('T2');

subplot(2,3,3);
imshow(PD, []);
title('Proton Density');

% 각 tissue mask 표시 (추가적으로 확인용)
subplot(2,3,4);
imshow(skull_mask, []);
title('Skull');

subplot(2,3,5);
imshow(grey_mask, []);
title('Grey Matter');

subplot(2,3,6);
imshow(white_mask | vent_mask, []); % white와 ventricle 포함
title('White Matter & Ventricle');

%% --- 함수 정의 ---
% ellipse_mask: (X,Y) 좌표에서 중심 (x0,y0), 반지름 a, b, 회전각 phi(도)를 적용한 타원 마스크 생성
function mask = ellipse_mask(X, Y, x0, y0, a, b, phi)
    % 회전각을 라디안으로 변환
    phi = deg2rad(phi);
    % 회전된 좌표 계산
    Xr = (X - x0) * cos(phi) + (Y - y0) * sin(phi);
    Yr = -(X - x0) * sin(phi) + (Y - y0) * cos(phi);
    % 타원 방정식 적용
    mask = ((Xr ./ a).^2 + (Yr ./ b).^2) <= 1;
end
