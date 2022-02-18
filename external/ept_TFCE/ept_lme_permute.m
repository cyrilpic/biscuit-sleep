function perm_table = ept_lme_permute(table, factors, participant_field)
% permutes the input table according to factors indicated and type of 
% randomization permissible (between or within randomisation)

% inputs
% ''''''
%
% 'table' should be the full data table used for observed analysis
% this table must include an identifier for each participant "participant_id"
%
% 'factors' is a struct containing two fields, 'name' and 'flag_within'
%    'name' is a cell array containing the names of the fields to be randomised
%    'flag_within' is a logical indicating which (if any) of the factors is a 
%     to be permuted within subject
% for example:
% factors = struct(...
%     'name', {'stroke_id', 'spindle_type'}, ... % covariate names
%     'flag_within', {0, 1}); % flag_within
%
% Original work from: https://github.com/Mensen/ept_TFCE-matlab
% Changelog (by @cyrilpic):
% - Make the participant field a parameter
% - Permutate between factors with the help of a group summary table
%   => accomodates all types of data

if nargin < 2
    error("missing factor structure")
end


between_variables_names = {factors(find(~[factors.flag_within])).name};
within_variables_names = {factors(find([factors.flag_within])).name};

summary = groupsummary(table, participant_field, @first, between_variables_names);
summary.GroupCount = [];
summary.Properties.VariableNames = [participant_field, between_variables_names(:)'];

% copy the table
perm_table = table;
perm_table(:, between_variables_names) = []; % remove between columns to be permutated;


% between factor randomisation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% find the unique group values
if ~isempty(between_variables_names)
    for b = 1 : length(between_variables_names)
        summary(:, between_variables_names{b}) = summary(randperm(height(summary)), between_variables_names{b});
    end
    perm_table = join(perm_table, summary);
end

% within factor randomisation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~


if ~isempty(within_variables_names)
    for v = 1 : length(within_variables_names)
        for n = 1 : length(summary.(participant_field))
            
            % get the relevant rows
            if ischar(summary.(participant_field)(n)) || iscellstr(summary.(participant_field)(n))
                rows = strcmp(table.(participant_field), summary.(participant_field)(n));
            else
                rows = (table.(participant_field) == summary.(participant_field)(n));
            end
            
            randomised_trials = randsample(table.(within_variables_names{v})(rows), ...
                length(table.(within_variables_names{v})(rows)), 0);
            
            % assign the new value
            perm_table.(within_variables_names{v})(rows) = randomised_trials;
            
        end
    end
end

end

% Sub function for group summary
function f = first(x)

if nargin == 0 || isempty(x)
    f = NaN;
else
    f = x(1);
end
end
