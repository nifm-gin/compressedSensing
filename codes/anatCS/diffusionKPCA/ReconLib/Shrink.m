

function y = Shrink( x, r )

a = abs(x);
m = ( a > 0 );
f = x ./ ( a + eps );

y = f .* m .* max( a - r, 0 );