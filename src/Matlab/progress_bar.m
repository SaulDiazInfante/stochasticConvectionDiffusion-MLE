function progress_bar(i, n, width, msg, f_bar)
    % progress_bar(i, n, width, msg)
    % i:    current iteration
    % n:    total iterations
    % width: width of the progress bar (default: 50)
    % msg:  optional user message to display

    if nargin < 3 || isempty(width)
        width = 50;
    end
    if nargin < 4
        msg = '';
    end

    percent = i / n;
    completed = floor(percent * width);
    remaining = width - completed;
    bar = [repmat('=', 1, completed), '>', repmat(' ', 1, max(0, remaining - 1))];
    rest = n - i;

    %fprintf(
    %    '\r[%s] %3.0f%% | Rem. : %d of %d | %s',
    %    bar,
    %    percent * 100,
    %    rest,
    %    n,
    %    msg
    %);
    %fflush(stdout);
    waitbar(percent, f_bar,...
        sprintf('\r[%s] %3.0f%% | Rem. : %d of %d | %s',...
            bar, ...
            percent * 100, ...
            rest, ...
            n, ...
            msg ...
        ) ...
    )
    %if i == n
    %    fprintf('\n');  % Newline at 100%
    %end
end
