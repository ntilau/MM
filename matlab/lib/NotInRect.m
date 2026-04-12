function bool = NotInRect( xmean, ymean, a, b, xo, yo)
% Test whether a point lies outside all rectangles in a set.
% ------------------------------------------------------------------------
% [IN]
%   xmean, ymean - query-point coordinates
%   a, b         - rectangle widths and heights
%   xo, yo       - rectangle centre offsets
%
% [OUT]
%   bool         - 1 if the point is outside every rectangle, 0 otherwise
%
set=0;
[x_vector, x_index]= sort(xo);
for i=1:length(x_index)-1
    idx=x_index(i);
    xmax(i) = xo(idx) + a(idx)/2;
    ymax(i) = yo(idx) + b(idx)/2;

    idx=x_index(i+1);
    xmin(i+1) = xo(idx) - a(idx)/2;
    ymin(i+1) = yo(idx) - b(idx)/2;

    xm=(xmax(i)+xmin(i+1))/2;
    ym=(ymax(i)+ymin(i+1))/2;

    if(xm == xmean && ym == ymean)
        set=1;
    end
end

if(set)
    bool = 1;
else
    bool = 0;
end
