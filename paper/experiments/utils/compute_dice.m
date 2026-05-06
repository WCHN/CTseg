function d = compute_dice(a, b)
% Dice coefficient between two binary volumes.
    sa = sum(a(:));
    sb = sum(b(:));
    if sa == 0 && sb == 0
        d = 1;  % both empty: perfect agreement
        return;
    end
    d = 2 * sum(a(:) & b(:)) / (sa + sb);
end
