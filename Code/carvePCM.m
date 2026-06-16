function [y_res, v_res] = carveMe2Modes_v2(model, M, eps, ...
    min_gr_photo, min_gr_hetero, rxn_scores, ...
    media, useGurobi, presolve, FeasibilityTol, saveName, solve)
    %description: MILP of the famous carveMe algorithm. Maximizes
    % rxn_scores while keeping biomass production at least at min_gr, for
    % both exchange media given via the media parameter. This is achieved
    % by duplicating the models and adding another binary variable to link
    % the two.
    % Obtains a flux-consistent network. 
    % input:
    %   model: Cobra model to be carved. Needs to contain c field to
    %       identify biomass
    %   M: upper/lower bound (1000)
    %   eps: min flux (1e-4)
    %   min_gr: min expected growth for photo mode
    %   min_gr: min expected growth for hetero mode
    %   rxn_scores: array of scores of length of rxns, as in carveME
    %   media: a table with columns 'rxn' (containing {'rxnId'}), 'photo',
    %   'hetero', indicating upper bound for the specified import rxns
    %   useGurobi: bool; if true, use gurobi as a solver
    %   presolve: bool; if true, first solve an LP to reduce num of vars
    %   FeasibilityTol: double, sets the gurobi FeasibilityTol param
    %   saveName: name of .mps file to write for solving the model on hpc
    %       (empty; do not write)
    %   solve: bool; if false, will only write the model and not solve it
    %
    % output:
    %   fluxes: table with fluxes of the obtained solution
    %   y_res: y values: 1 indicates that a reaction is kept. Consists of
    %       y_f, y_b (first forward, then backward for all rxns)
    %   v_res: flux values of the obtained solution

    if nargin < 11
        saveName = '';
    end
    if nargin < 12
        solve = true;
    end

    % Extract stoichiometric matrix and problem dimensions
    S = full(model.S);  % Convert sparse to full if needed
    [Nm, N] = size(S);
    bio_ix = find(model.c, 1);

    % Define bounds
    v_lb_photo = [model.lb];
    v_lb_hetero = [model.lb];
    v_lb_photo(bio_ix) = min_gr_photo;
    v_lb_hetero(bio_ix) = min_gr_hetero;

    v_ub_photo = [model.ub];
    v_ub_hetero = [model.ub];
    for rxnIx = 1:length(media.rxn)
        currRxnIx = strcmp(media.rxn{rxnIx}, model.rxns);
        if any(currRxnIx)
            v_ub_photo(currRxnIx) = media.photo(rxnIx);
            v_ub_hetero(currRxnIx) = media.hetero(rxnIx);
        end
    end

    % display different bounds between growth modes
    diffIdx = find(v_ub_photo ~= v_ub_hetero);
    if ~isempty(diffIdx)
        disp('Differing upper bounds between photo and hetero:')
    end
    for i = 1:length(diffIdx)
        fprintf([model.rxns{diffIdx(i)}, '\t', ...
                num2str(v_ub_photo(diffIdx(i))), '\t', ...
                num2str(v_ub_hetero(diffIdx(i))), '\n']);
    end
    diffIdx = find(v_lb_photo ~= v_lb_hetero);
    if ~isempty(diffIdx)
        disp('Differing lower bounds between photo and hetero:')
    end
    for i = 1:length(diffIdx)
        fprintf([model.rxns{diffIdx(i)}, '\t', ...
            num2str(v_lb_photo(diffIdx(i))), '\t', ...
            num2str(v_lb_hetero(diffIdx(i))), '\n']);
    end

    % Set up variables and objective function
    % variables: y, y_fp, y_rp, yfh, yrh, v_p, v_h
    f = [-rxn_scores; zeros(4 * N, 1); zeros(2 * N, 1)]; % c * y
    intcon = 1:(5 * N);  % Integer variables
    lb = [zeros(5 * N, 1); v_lb_photo(:); v_lb_hetero(:)];
    ub = [ones(5 * N, 1); v_ub_photo(:); v_ub_hetero(:)];
    
    % Steady state constraints
    % S*v = 0 for steady state for each v_p and v_h
    Aeq = [sparse(Nm, (5 * N)), S, sparse(Nm, N); 
           sparse(Nm, (5 * N)), sparse(Nm, N), S];
    beq = zeros(2 * Nm, 1);

    % Construct MILP constraints
    eyeN = speye(N);
    Aineq1 = [sparse(N, N), eyeN * eps, eyeN * -M, sparse(N, 2 *N), -eyeN, sparse(N, N)];
    Aineq2 = [sparse(N, 3 * N), eyeN * eps, eyeN * -M, sparse(N, N), -eyeN];
    Aineq3 = [sparse(N, N), eyeN * -M, eyeN * eps, sparse(N, 2 *N), eyeN, sparse(N, N)];
    Aineq4 = [sparse(N, 3 * N), eyeN * -M, eyeN * eps, sparse(N, N), eyeN];
    Aineq5 = [sparse(N, N), eyeN, eyeN, sparse(N, 4 * N)];
    Aineq6 = [sparse(N, 3 * N), eyeN, eyeN, sparse(N, 2 * N)];
    Aineq7 = [eyeN * -2, eyeN, eyeN, eyeN, eyeN, sparse(N, 2 * N)];
    bineq1 = zeros(N, 1);  % rhs must be dense for gurobi
    bineq2 = zeros(N, 1);
    bineq3 = zeros(N, 1);
    bineq4 = zeros(N, 1);
    bineq5 = ones(N, 1);
    bineq6 = ones(N, 1);
    bineq7 = zeros(N, 1);

    % Combine inequality constraints
    Aineq = [Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6; Aineq7];
    bineq = [bineq1; bineq2; bineq3; bineq4; bineq5; bineq6; bineq7];

    % Solve MILP problem using MATLAB's solver
    % much slower the first time I tried it (>2hrs vs. 5 seconds)
    if ~useGurobi
        options = optimoptions('intlinprog', 'Display', 'off');
        [x, fval, exitflag, output] = intlinprog(f, intcon, Aineq, ...
            bineq, Aeq, beq, lb, ub, options);
        % Check optimal status
        if exitflag == 1
            v_res = x((5 * N + 1):end);
            y_res = x(1:5 * N);
            
            return
        else
            disp('Model not optimal. Status:');
            disp(exitflag);
            disp('fval:');
            disp(fval);
            disp('output stored in y_res, x in v_res');
            v_res = x;
            y_res = output;
            
            return
        end
    else
        % instead, use gurobi
        gModel = struct();
        gModel.A = sparse([Aineq;Aeq]);
        gModel.obj = f;
        gModel.rhs = [bineq;beq];
        gModel.sense = [repmat('<', 7*N, 1);
                        repmat('=', 2*Nm, 1)];
        gModel.lb = lb;
        gModel.ub = ub;
        gModel.modelsense = 'min';
        
        params = struct();
        params.FeasibilityTol = FeasibilityTol;
        params.IntFeasTol = FeasibilityTol;
        params.OptimalityTol = FeasibilityTol;

        if presolve  % not sure this makes a lot of sense
            changed = true;
            while changed
                changed = false;
                gModel.vtype = [repmat('c', 5*N, 1);
                                repmat('c', 2*N, 1)];
                result = gurobi(gModel, params);
                y_res = result.x(1:(5 * N));
                to_1 = false;
                % we can only fix variables if y is 1
                for i = 1:length(y_res)
                    if lb(i) == 1
                        continue
                    end
                    if y_res(i) + 1e-7 > 1
                        lb(i) = 1;
                        to_1 = to_1 + 1;
                        changed = true;
                    end
                end
            end

            disp(['Presolve resulted in ' num2str(to_1) ...
                ' variables fixed to 1']);
        end
        gModel.vtype = [repmat('B', 5*N, 1);
                        repmat('c', 2*N, 1)];
        if ~isempty(saveName)
            gurobi_write(gModel, [saveName '.mps']);
        end
        if ~solve
            y_res = [];
            v_res = [];
            return;
        end
        result = gurobi(gModel, params);

        if strcmp(result.status, 'OPTIMAL')
            y_res = result.x(1:(5 * N));
            v_res = result.x((5 * N + 1):end);
            disp(['Optimal objective: ', num2str(result.objval)]);
        else
            disp(['Optimization returned status: ', ...
                num2str(result.status)]);
            disp('result struct is stored in y_res');
            y_res = result;
            v_res = [];
        end
    end
end