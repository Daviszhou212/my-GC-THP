%%% DFT codebook genetation for UPA
function [codebook, CB_size] = dftCodebook(M, N, verAng, horAng)
    % M and N are # of vertical and horizontal elements respectively on each panel
    % verAng and horAng are ver/hor angles between beams and panels' normal line in [degree]

    temp = 2 * (1:M) / M;
    temp(temp > 1) = temp(temp > 1) - 2;
    verCos = temp(abs(temp) <= cos(pi * (90 - verAng) / 180));

    temp = 2 * (1:N) / N;
    temp(temp > 1) = temp(temp > 1) - 2;
    horCos = temp(abs(temp) <= cos(pi * (90 - horAng) / 180));

    CB_size = length(verCos) * length(horCos);
    codebook = zeros(M * N, CB_size);
    index = 1;

    for x = verCos

        for y = horCos
            codebook(:, index) = kron(exp(1j * y * pi * (0:N - 1)).', exp(1j * pi * x * (0:M - 1)).');
            index = index + 1;
        end

    end

    codebook = codebook ./ sqrt(M * N);
end