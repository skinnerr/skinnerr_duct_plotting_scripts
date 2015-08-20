function [ position ] = Centered_Figure_Position( figure_width, figure_height )
%%%
% Calculates the position array necessary to center a figure of given height and width on
% the root display.
%%%

    fw = figure_width;
    fh = figure_height;
    screen_size = get(groot,'ScreenSize');
    sw = screen_size(3);
    sh = screen_size(4);

    assert(fw > 0, 'Figure width must be greater than 0.');
    assert(fh > 0, 'Figure height must be greater than 0.');
    
    % If figure would be larger than the display, resize it to fit with some padding.
    if fw > sw
        fw = ceil(0.9 * sw);
    end
    if fh > sh
        fh = ceil(0.9 * sh);
    end

    % Calculate left and bottom offsets that result in a centered figure.
    fl = (sw - fw) / 2;
    fb = (sh - fh) / 2;
    
    position = [fl, fb, fw, fh];

end

