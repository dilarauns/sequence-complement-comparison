

function found = checkUp(temp1, temp2, minHarf, N)
% 



%     while 1
%         
%             disp('%%%%%%%%%%%% check updan gelen sonuç %%%%%%%%%%%%%%%') ;
%             fprintf('Min harf : %d \n', minHarf) ; 
%             
    for k=1:(N - minHarf)
       for a = 1:length(temp1) - (minHarf - 1)
           
            disp('%%%%%%%%%%%% check updan gelen sonuç %%%%%%%%%%%%%%%') ;
            found = strfind(temp2,temp1(a : a + (minHarf - 1))) ; 
            fprintf('Temp 1 : %s \n', temp1) ; 
            %disp(temp1);
            fprintf('Found : %s \n',temp2(found : found + (minHarf - 1)));
            fprintf('Temp 2 : %s \n', temp2) ; 
            fprintf('Found Index : %d \n', found) ;
            disp('%%%%%%%%%%%% check updan gelen sonuç %%%%%%%%%%%%%%%') ;
            % İLK BULUŞTAN SONRA KONTROL
            
            if (isempty(found) ==0)
                
                minHarf = minHarf + 1 ; 
            
            else 
                break ; 
            end
       end % for 
     end % for k

%         
%     end % end while   

    
end % function end