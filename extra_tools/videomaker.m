%Place in directory of image sequence and run
Title_vid='0deg.avi'; % select video title
image_names = dir('*.png'); % pattern to match filenames.
%optional for title
title_c=0;
if title_c==0
    tmp_fig=figure();
    mkdir temp_img
    cd temp_img 
    movie_obj = VideoWriter(Title_vid);
    open(movie_obj);
    for K = 1 : length(image_names)
      this_image = imread(['../' image_names(K).name]);
      imshow(this_image)
      title(image_names(K).name, 'interpreter', 'none');
      saveas(gcf,[sprintf( '%04d', K ) '.png'])
      image_names_tmp = dir('*.png');
      this_image_tmp = imread(image_names_tmp(K).name); 
      for i=1:2
        writeVideo(movie_obj, this_image_tmp);
      end
    end
else   
    movie_obj = VideoWriter(Title_vid);
    open(movie_obj);
    for K = 1 : length(image_names)
      this_image = imread(image_names(K).name);
      for i=1:2
        writeVideo(movie_obj, this_image);
      end
    end
end
close(movie_obj);
%implay(Title_vid)


