clc ; 
clear ;

decision = input("Dizinizi girmek için 1, Rastgele dizi oluşturmak için 2 : " );

while (1)
    if (decision == 1) 
            
        dna = upper(input("Dizinizi giriniz : ", "s")) ;
        break ;
    elseif (decision == 2)
    
        base = ['A', 'C', 'G', 'T'] ; 
    
        dnaL = input("DNA UZUNLUĞU GIRINIZ : ") ;
        dnaIndexes = randi([1, length(base)], 1, dnaL) ; 
        dna = base(dnaIndexes) ; 
        break ; 
    else
        disp('Geçersiz değer girdiniz.') ;   
        decision = input("Dizinizi girmek için 1, Rastgele dizi oluşturmak için 2 : " );
    end % if end 
end % while end

%% section deneme
N = input("Sekans sayısı giriniz : ") ; 
NL = input("Sekans uzunluğu giriniz :") ;

%dna =  "ATCGATACGCTAGCATGCAGCATCAGCATCAGATCGATCACGACTACCGACTACAGCTACCGATCGCATCACGATCAGATCAGCATCAGCTACAGCTAACTACGCATCAGCATCGACTACACACACGATCGATA" ; 

randI = randperm(length(dna) - (NL - 1)) ;
sequenceArray = strings(1, N) ;
complementArray = strings(1,N) ;
 
for i=1:N 
    
    sequenceArray(i) = dna(randI(i): randI(i) + (NL-1) ) ; 
    temp = cell2mat(sequenceArray(i)) ;
    for z = 1:NL
         
        %disp(temp);
        if temp(z) == 'A' 
            temp(z) = 'T' ;
        elseif temp(z) == 'T'
            temp(z) = 'A';
        elseif temp(z) == 'G' 
            temp(z) = 'C';
        elseif temp(z) == 'C' 
            temp(z) = 'G' ;
        end
    end
    
    %disp(flip(temp));
    complementArray(i) = flip(temp) ; 
end % end for 

disp("Sequence Array Complement Array" );
disp([ sequenceArray' , complementArray']);




indel = input('Indel puanını giriniz : ');
mismatch = input('Mismatch puanını giriniz : ');
match = input('Match puanını giriniz : ') ;
score = input('En küçük skor değerini giriniz : ') ; 



overlapScoreMatrix = ones(N, N) * -99 ; 
matchSequenceMatrix = strings(N, N);

%cellMatchSequenceMatrix = cell(matchSequenceMatrix); 
file = fopen('layout.txt', 'w');
fprintf(file,'SEQUENCE - SEQUENCE KARŞILAŞTIRMASI\n');
count = 0 ; 
for a = 1:length(sequenceArray) 
    temp1 = cell2mat(sequenceArray(a)) ;
    for b = a+1 :length(sequenceArray)
        temp2 = cell2mat(sequenceArray(b)) ;
        
        fprintf("Temp1 : %s \n", temp1) ; 
        fprintf("Temp2 : %s \n", temp2) ;
        readReadMatrix = zeros(NL + 1, NL + 1) ; 
        for i= 2 : length(temp1) + 1
            for j = 2 : length(temp1) + 1
    
                if (temp1(i - 1) == temp2( j - 1 ))
                    tempScore = match ; 
                else
                    tempScore = mismatch * -1 ;
                end

                readReadMatrix(i,j) = max([readReadMatrix(i-1, j) - indel, readReadMatrix(i-1, j-1) + tempScore, readReadMatrix(i, j-1) - indel ]);
                
                

                
                
            end % for k 
        
                
            
            disp(readReadMatrix) ;
            [bestOverlap,column] = max(max(readReadMatrix)) ;
            [~,row] = max(readReadMatrix(:,column));
            %fprintf('column : %d, row : %d \n', column, row) ;           
        
        end % for i 
        

%         fprintf('Best overlap %d \n', bestOverlap) ;
        if bestOverlap >= score 
            overlapScoreMatrix(a,b) = bestOverlap ; 
            

            % BURANIN AMACI ARANACAK STRINGI VE ARANACAK INDISIN BAŞLANGIÇ
            % NOKTASINI BULMAK.
            temp1Str = "" ;
            temp2Str = "" ;
            

            i = row ; 
            j = column ; 
            while ( i ~= 1  &&  j~=1 )
                    %fprintf('i : %d, j : %d\n', i, j) ;
                    
                    
                    if readReadMatrix(i, j-1) == readReadMatrix(i,j) + indel
                        % SOL
                        temp1Str = append(temp1Str, '-') ; 
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        %fprintf('TEMP1 INDEL  \n');
                        if j ~= 1
                            j = j-1 ; 
                        end
                   
                    elseif readReadMatrix(i-1, j) == readReadMatrix(i,j) + indel
                        % AŞAĞI YUKARI
                        temp1Str = append(temp1Str, temp1(i-1)) ; 
                        temp2Str = append(temp2Str, '-') ;  
                        %fprintf('TEMP2 INDEL  \n');
                        
                         if i ~= 1
                            i = i-1 ;
                            
                         end
                        
                    
                    elseif readReadMatrix(i-1, j-1) == readReadMatrix(i,j) - match
                        % ÇAPRAZ
                        temp1Str = append(temp1Str,temp1(i-1)) ;
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        %fprintf('MATCH OLMUŞ MU i,j : %d,%d temp1 %s ,temp1Check %s , kontrol : %s\n',i,j,temp1 ,temp1(i-1), temp1Str );
                        %fprintf('MATCH OLMUŞ MU i,j : %d,%d temp2 %s ,temp2Check %s , kontrol : %s\n',i,j,temp2 ,temp2(j-1), temp2Str );
                        
                        if i ~= 1
                            i = i-1 ;
                        end

                        if j ~= 1
                            j = j-1 ; 
                        end
    
                    
                    elseif readReadMatrix(i-1, j-1) == readReadMatrix(i,j) + mismatch
                        % MISMATCH
                        temp1Str = append(temp1Str,temp1(i-1)) ; 
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        %fprintf('MISMATCH OLMUŞ MU i,j : %d,%d temp1 %s ,temp1Check %s , kontrol : %s\n',i,j,temp1 ,temp1(i-1), temp1Str );
                        %fprintf('MISMATCH OLMUŞ MU i,j : %d,%d temp2 %s ,temp2Check %s , kontrol : %s\n',i,j,temp2 ,temp2(j-1), temp2Str );
                        if i ~= 1
                            i = i-1 ;
                        end

                        if j ~= 1
                            j = j-1 ; 
                        end
                     
                    else
                        % BAŞLANGICA DÖNMÜŞ
                  
                        %cellMatchSequenceMatrix(a,b) = {temp1Str , temp2Str} ; 
                        
                        if i ~= 1
                            i = i-1 ;
                        end

                        if j ~= 1
                            j = j-1 ; 
                        end

                        if i == 1 && j == 1
                            break ; 
                        end


            
            
                    end % end first if 


            end % end while
          
            temp1Str = flip(cell2mat(temp1Str)) ;
            temp2Str = flip(cell2mat(temp2Str)) ;
            fprintf('temp1Str : %s \n', temp1Str) ; 
            fprintf('temp2Str : %s \n', temp2Str) ; 

            if (i ~= 2)
                
                unMatchedStr = temp1(1:i-1);
                unMatchedStr = append(unMatchedStr, temp1Str) ; 
                
                if (row ~= NL + 1)
                    
                    unMatchedStr = append(unMatchedStr, temp1(row:length(temp1)));
        
                else 
                end
            
                %fprintf('Hizalanmış   : %s\n', unMatchedStr);
            elseif ( i == 2 )

                unMatchedStr = temp1(1) ;

                unMatchedStr = append(unMatchedStr, temp1Str) ;
                
                unMatchedStr = append(unMatchedStr, temp1(row:length(temp1)));
                %fprintf('Hizalanmış   : %s\n', unMatchedStr);
            end %% temp1 str kontrol if i 
            


            %temp2str
            if (j ~= 2)
                
                unMatchedStr2 = temp2(1:j-1);
                unMatchedStr2 = append(unMatchedStr2, temp2Str) ; 
                
                if (column ~= NL + 1)
                    
                    unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
        
                
                

                else 
                end
                %fprintf('Hizalanmış 2 : %s\n', unMatchedStr2);

            elseif ( j == 2 )
                unMatchedStr2 = temp2(1) ;

                unMatchedStr2 = append(unMatchedStr2, temp2Str) ;
                
                unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
                %fprintf('Hizalanmış 2 : %s\n', unMatchedStr2);
            end %% temp1 str kontrol if i 

            
            if (i > j)
                
                fprintf(file,'Hizalanmış 1 : %s\n', unMatchedStr);
                fprintf(file,'Hizalanmış 2 : ');
                fprintf(file,repmat(' ',1,i-j));
                fprintf(file,'%s\n',unMatchedStr2);
                fprintf(file,'------------------------------------------------\n');
            
            elseif(j>i)
                
                fprintf(file,'Hizalanmış 1 : ');
                fprintf(file,repmat(' ',1,j-i));
                fprintf(file,'%s\n',unMatchedStr);
                fprintf(file,'Hizalanmış 2 : %s\n',unMatchedStr2);
                fprintf(file,'------------------------------------------------\n');
                
                

             elseif(i==j)
                 fprintf(file,'Hizalanmış 1 : %s\n', unMatchedStr);
                 fprintf(file,'Hizalanmış 2 : %s\n', unMatchedStr2);
                 fprintf(file,'------------------------------------------------\n');
             end


        disp('------------------------------------------------');
            
        end % end bestoverlap if

 
    end % for b 
    
end % for a 
% bestOverlap = readReadMatrix(length(sequenceArray), length(sequenceArray)) ;






% READ COMPLEMENT 
fprintf(file, 'SEQUENCE COMPLEMENT KARŞILAŞTIRMASI\n');
for x = 1 : length(complementArray)
    temp1 = cell2mat(complementArray(x)) ;
    for y = x+1 : length(complementArray)
        temp2 = cell2mat(sequenceArray(y)) ;
%         fprintf("Seq : %s \n", seq) ; 
%         fprintf("Comp : %s \n", comp) ;
        
        readComplementMatrix = zeros(NL + 1, NL + 1) ; 
        
        for k = 2 : length(temp1) + 1
            for d = 2 : length(temp2) + 1
                if (temp1(k - 1) == temp2( d - 1 ))
                    score2 = match ; 
                else
                    score2 = mismatch * -1 ; 
                end
                readComplementMatrix(k,d) = max([readComplementMatrix(k-1, d) - indel, readComplementMatrix(k-1, d-1) + score2, readComplementMatrix(k, d-1) - indel ]);
%                 disp(readComplementMatrix) ;
            end % for d 
            [bestCompOverlap,column] = max(max(readComplementMatrix)) ;
            [~,row] = max(readComplementMatrix(:,column));
            %fprintf('column : %d, row : %d \n', column, row) ;
        end % for k
%         disp('------------------------------------------------');
%      




        if bestCompOverlap >= score 
            overlapScoreMatrix(y,x) = bestCompOverlap ; 
            
    
            % BURANIN AMACI ARANACAK STRINGI VE ARANACAK INDISIN BAŞLANGIÇ
            % NOKTASINI BULMAK.
            temp1Str = "" ;
            temp2Str = "" ;
            boolTemp = 0 ;
    
            i = row ; 
            j = column ; 
            while ( i ~= 1  &&  j~=1 )
                    %fprintf('i : %d, j : %d\n', i, j) ;
                    
                    
                    if readComplementMatrix(i, j-1) == readComplementMatrix(i,j) + indel
                        % SOL
                        temp1Str = append(temp1Str, '-') ; 
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        %fprintf('TEMP1 INDEL  \n');
                        if j ~= 1
                            j = j-1 ; 
                        end
                   
                    elseif readComplementMatrix(i-1, j) == readComplementMatrix(i,j) + indel
                        % AŞAĞI YUKARI
                        temp1Str = append(temp1Str, temp1(i-1)) ; 
                        temp2Str = append(temp2Str, '-') ;  
                        %fprintf('TEMP2 INDEL  \n');
                        
                         if i ~= 1
                            i = i-1 ;
                            
                         end
                        
                    
                    elseif readComplementMatrix(i-1, j-1) == readComplementMatrix(i,j) - match
                        % ÇAPRAZ
                        temp1Str = append(temp1Str,temp1(i-1)) ;
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        %fprintf('MATCH OLMUŞ MU i,j : %d,%d temp1 %s ,temp1Check %s , kontrol : %s\n',i,j,temp1 ,temp1(i-1), temp1Str );
                        %fprintf('MATCH OLMUŞ MU i,j : %d,%d temp2 %s ,temp2Check %s , kontrol : %s\n',i,j,temp2 ,temp2(j-1), temp2Str );
                        
                        if i ~= 1
                            i = i-1 ;
                        end
    
                        if j ~= 1
                            j = j-1 ; 
                        end
    
                    
                    elseif readComplementMatrix(i-1, j-1) == readComplementMatrix(i,j) + mismatch
                        % MISMATCH
                        temp1Str = append(temp1Str,temp1(i-1)) ; 
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        %fprintf('MISMATCH OLMUŞ MU i,j : %d,%d temp1 %s ,temp1Check %s , kontrol : %s\n',i,j,temp1 ,temp1(i-1), temp1Str );
                        %fprintf('MISMATCH OLMUŞ MU i,j : %d,%d temp2 %s ,temp2Check %s , kontrol : %s\n',i,j,temp2 ,temp2(j-1), temp2Str );
                        if i ~= 1
                            i = i-1 ;
                        end
    
                        if j ~= 1
                            j = j-1 ; 
                        end
                     
                    else
                        % BAŞLANGICA DÖNMÜŞ
                  
                        %cellMatchSequenceMatrix(a,b) = {temp1Str , temp2Str} ; 
                        
                        if i ~= 1
                            i = i-1 ;
                        end
    
                        if j ~= 1
                            j = j-1 ; 
                        end
    
                        if i == 1 && j == 1
                            break ; 
                        end
    
    
            
            
                    end % end first if 
    
    
            end % end while
          
            temp1Str = flip(cell2mat(temp1Str)) ;
            temp2Str = flip(cell2mat(temp2Str)) ;
            fprintf('temp1Str : %s \n', temp1Str) ; 
            fprintf('temp2Str : %s \n', temp2Str) ; 
    
            if (i ~= 2)
                
                unMatchedStr = temp1(1:i-1);
                unMatchedStr = append(unMatchedStr, temp1Str) ; 
                
                if (row ~= NL + 1)
                    
                    unMatchedStr = append(unMatchedStr, temp1(row:length(temp1)));
        
                else 
                end
    
                %fprintf('Hizalanmış   : %s\n', unMatchedStr);
            elseif ( i == 2 )
    
                unMatchedStr = temp1(1) ;
    
                unMatchedStr = append(unMatchedStr, temp1Str) ;
                
                unMatchedStr = append(unMatchedStr, temp1(row:length(temp1)));
                %fprintf('Hizalanmış   : %s\n', unMatchedStr);
            end %% temp1 str kontrol if i 
            
    
    
            %temp2str
            if (j ~= 2)
                
                unMatchedStr2 = temp2(1:j-1);
                unMatchedStr2 = append(unMatchedStr2, temp2Str) ; 
                
                if (column ~= NL + 1)
                    
                    unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
        
                
                
    
                else 
                end
                %fprintf('Hizalanmış 2 : %s\n', unMatchedStr2);
    
            elseif ( j == 2 )
                unMatchedStr2 = temp2(1) ;
    
                unMatchedStr2 = append(unMatchedStr2, temp2Str) ;
                
                unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
                %fprintf('Hizalanmış 2 : %s\n', unMatchedStr2);
            end %% temp1 str kontrol if i 
    
            
            if (i > j)
                
                fprintf(file,'Hizalanmış 1 : %s\n', unMatchedStr);
                fprintf(file,'Hizalanmış 2 : ');
                fprintf(file,repmat(' ',1,i-j));
                fprintf(file,'%s\n',unMatchedStr2);
                fprintf(file,'------------------------------------------------\n');
            
            elseif(j>i)
                
                fprintf(file,'Hizalanmış 1 : ');
                fprintf(file,repmat(' ',1,j-i));
                fprintf(file,'%s\n',unMatchedStr);
                fprintf(file,'Hizalanmış 2 : %s\n',unMatchedStr2);
                fprintf(file,'------------------------------------------------\n');
                
                

             elseif(i==j)

                 fprintf(file,'Hizalanmış 1 : %s\n', unMatchedStr);
                 fprintf(file,'Hizalanmış 2 : %s\n',unMatchedStr2);
                 fprintf(file,'------------------------------------------------\n');
                 
             end


        disp('------------------------------------------------');
        
        end % end bestoverlap if


    end % for y 
end % for x 

     
fclose(file);

















