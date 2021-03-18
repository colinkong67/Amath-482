[images, labels] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
images = double(reshape(images, size(images,1)*size(images,2), []).');
labels = double(labels);

[testdata, testgnd] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');
testdata = double(reshape(testdata, size(testdata,1)*size(testdata,2), []).');
testgnd = double(testgnd);
[row,col,num] = size(images);
reshape_images = zeros(row*col,num);
for i=1:num
    reshape_images(:,i) = reshape(images(:,:,i),row*col,1);
end

[U,S,V] = svd([testdata testgnd],'econ');
for k = 1:10
    subplot(2,5,k)
    ut1 = reshape(U(:,k),100,100);
    ut2 = rescale(ut1);
    imshow(ut2)
end




