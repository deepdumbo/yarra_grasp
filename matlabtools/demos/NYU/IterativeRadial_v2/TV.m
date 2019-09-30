% ## Helper function for calculating the TV and its gradient for a given
% ## 2D complex image. TV is calculated using only the real part, and using
% ## first-order finite differences.

function [f,g]=TV(img)
    % Assuming that the image is squared
    br=size(img,1);

    % Initialize the output variables (f is double, and g complex matrix of
    % size br x br
    f=0;
    g=complex(zeros(br),0);
    
    % Loop over image and simultanously calculate the penalty and the
    % gradient. The edge pixels are excluded here for simplicity (a better
    % solution would be to "connect" the left-rigt and top-bottom edges)
    for ix=2:br-1
        for iy=2:br-1
            f=f+abs(real(img(ix,iy))-real(img(ix-1,iy)));
            f=f+abs(real(img(ix,iy))-real(img(ix,iy-1)));

            gradval=0;        

            gradval=gradval+sign( real(img(ix,iy))-real(img(ix-1,iy)) );
            gradval=gradval-sign( real(img(ix+1,iy))-real(img(ix,iy)) );

            gradval=gradval+sign( real(img(ix,iy))-real(img(ix,iy-1)) );
            gradval=gradval-sign( real(img(ix,iy+1))-real(img(ix,iy)) );

            g(ix,iy)=gradval + 1i*0;
        end
    end
    
end

