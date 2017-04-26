function closestOdd = round2ClosestOddNum(x,roundingDirection)
if nargin<2
    roundingDirection = 'round';
end
switch roundingDirection
    case {'ceil','up'}
        closestOdd = 2*ceil((x-1)/2)+1;
    case {'floor','down'}
        closestOdd = 2*floor((x-1)/2)+1;
    case 'round'
        closestOdd = 2*round((x-1)/2)+1;
    otherwise
        error('Not supported');
end