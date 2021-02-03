x = LinRange(-100, 100, 100);
y = LinRange(-100, 100, 100);
X = [x, y];

arr = rand(100, 100) + im*rand(100, 100)
grad = gradient_2D(arr, X);

wps_v_2D(grad, [50, 50])