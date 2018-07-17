
function create_path_and_map ()

global POSES LM

POSES= [0,0,0; 10000,0,0];
% POSES= [0,0,pi/2; 25,25,0; 50,0,-pi/2; 25,-25,pi; 0,0,pi/2];
xCoord= linspace(0,30,2);
xCoord= [xCoord, xCoord];
yCoord= ones(1,4)*5;
yCoord(3:end)= yCoord(3:end)*(-1);

% posLM(1,:)= [-25,75];
% negLM(1,:)= [15,25];

% well spaced
% posLM(2,:)= [0,0];
% negLM(2,:)= [0,0];

% Bad spaced
% posLM(2,8:102)= linspace(-2.5,-2.5,102-7);
% negLM(2,8:102)= linspace(2.5,2.5,102-7);
LM= [xCoord;yCoord];












