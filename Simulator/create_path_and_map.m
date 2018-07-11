
function create_path_and_map ()

global POSES LM

% POSES= [0,0,0; 10000,0,0];
POSES= [0,0,pi/2; 25,25,0; 50,0,-pi/2; 25,-25,pi; 0,0,pi/2];
% posLM(1,:)= linspace(20,100,3);
% negLM(1,:)= linspace(20,100,3);
posLM(1,:)= [-25,75];
negLM(1,:)= [15,25];

% well spaced
posLM(2,:)= [0,0];
negLM(2,:)= [0,0];

% Bad spaced
% posLM(2,8:102)= linspace(-2.5,-2.5,102-7);
% negLM(2,8:102)= linspace(2.5,2.5,102-7);
LM= [posLM,negLM];












