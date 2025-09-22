% Copyright 2024
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [out]=eval_maxStress_COMSOL(Stress,mat)

stress_rot =[];

for kk = 1:size(Stress, 1)
    if Stress(kk, 3) > mat.Rotor.sigma_max
        stress_rot = [stress_rot; Stress(kk,1), Stress(kk,2)];
    end
end

[~, idx] = max(Stress(:, 3));

out.stressrot = stress_rot;

%out.x_over    = stress_rot(:,1);
%out.y_over    = stress_rot(:,2);
if ~isempty(stress_rot)
    out.x_over = stress_rot(:, 1);
    out.y_over = stress_rot(:, 2);
else
    out.x_over = [];
    out.y_over = [];
end

out.x_max     = Stress(idx, 1);
out.y_max     = Stress(idx, 2);
