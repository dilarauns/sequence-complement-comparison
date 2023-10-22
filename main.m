
clc ; 
clear ;
%kullanıcı karar aşaması

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


N = input("Sekans sayısı giriniz : ") ; 
NL = input("Sekans uzunluğu giriniz :") ;

%dna içerisinden rastgele parça alma 

randI = randperm(length(dna) - (NL - 1)) ;
sequenceArray = strings(1, N) ;
complementArray = strings(1,N) ;
 
for i=1:N 
    
    sequenceArray(i) = dna(randI(i): randI(i) + (NL-1) ) ; 
    temp = cell2mat(sequenceArray(i)) ;
    for z = 1:NL
         
      
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
    
    
    complementArray(i) = flip(temp) ; 
end % end for 

disp("Sequence Array Complement Array" );
disp([ sequenceArray' , complementArray']);



%kullanıcıdan skor için alınan bilgiler
indel = input('Indel puanını giriniz : ');
mismatch = input('Mismatch puanını giriniz : ');
match = input('Match puanını giriniz : ') ;
score = input('En küçük skor değerini giriniz : ') ; 



%overlap score matris ve eşleşme matrisi oluşturma 
overlapScoreMatrix = ones(N, N) * -99 ; 
matchSequenceMatrix = strings(N, N);

%cellMatchSequenceMatrix = cell(matchSequenceMatrix);

%layoutları dosyalama işlemi-open folder
file = fopen('layout.txt', 'w');
fprintf(file,'SEQUENCE - SEQUENCE KARŞILAŞTIRMASI\n');


fprintf('SEQUENCE - SEQUENCE KARŞILAŞTIRMASI\n');

%cell2mat -- sekansları string dizisine çevirme
%sırayla sekandların karşılaştırılması işlemi
for a = 1:length(sequenceArray) 
    temp1 = cell2mat(sequenceArray(a)) ;
    for b = a+1 :length(sequenceArray)
        temp2 = cell2mat(sequenceArray(b)) ;
        
        fprintf("Temp1 : %s \n",temp1) ; 
        fprintf("Temp2 : %s \n", temp2) ;
        readReadMatrix = zeros(NL + 1, NL + 1) ; 
       %cell2mat ile oluşturulan dizilerin karşılaştırılması
        for i= 2 : length(temp1) + 1
            for j = 2 : length(temp1) + 1
    
                if (temp1(i - 1) == temp2( j - 1 ))
                    tempScore = match ; 
                else
                    tempScore = mismatch * -1 ;
                end

                readReadMatrix(i,j) = max([readReadMatrix(i-1, j) - indel, readReadMatrix(i-1, j-1) + tempScore, readReadMatrix(i, j-1) - indel ]);
                
                

                
                
            end % for k 
        
                
            
            %matris maksimum değer tutucu
            [bestOverlap,column] = max(max(readReadMatrix)) ;
            [~,row] = max(readReadMatrix(:,column));
                      
        
        end % for i 
        

           %overlap matris işlemesi
        if bestOverlap >= score 
            overlapScoreMatrix(a,b) = bestOverlap ; 
            

            
            
            temp1Str = "" ;
            temp2Str = "" ;
          
           

            i = row ; 
            j = column ; 
            while ( i ~= 1  &&  j~=1 )
                   
                    % her sekans için eşleşme durumlarını matriste tersten giderek yazdırma
                    
                    if readReadMatrix(i, j-1) == readReadMatrix(i,j) + indel
                        % sağdan gelme durumu - indel for temp1
                        temp1Str = append(temp1Str, '-') ; 
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        
                        %while koşulları
                        if j ~= 1
                            j = j-1 ; 
                        end
                   
                    elseif readReadMatrix(i-1, j) == readReadMatrix(i,j) + indel
                        % yukardan gelme durumu-indelfor temp2
                        temp1Str = append(temp1Str, temp1(i-1)) ; 
                        temp2Str = append(temp2Str, '-') ;  
                        
                        %while koşulları
                         if i ~= 1
                            i = i-1 ;
                            
                         end
                        
                    
                    elseif readReadMatrix(i-1, j-1) == readReadMatrix(i,j) - match
                        % çaprazdan gelme durumu
                        temp1Str = append(temp1Str,temp1(i-1)) ;
                        temp2Str = append(temp2Str,temp2(j-1)) ;  

                        %kontrol satırları
                        %fprintf('MATCH OLMUŞ MU i,j : %d,%d temp1 %s ,temp1Check %s , kontrol : %s\n',i,j,temp1 ,temp1(i-1), temp1Str );
                        %fprintf('MATCH OLMUŞ MU i,j : %d,%d temp2 %s ,temp2Check %s , kontrol : %s\n',i,j,temp2 ,temp2(j-1), temp2Str );
                        
                        if i ~= 1
                            i = i-1 ;
                        end

                        if j ~= 1
                            j = j-1 ; 
                        end
    
                    
                    elseif readReadMatrix(i-1, j-1) == readReadMatrix(i,j) + mismatch
                        % mismatch durumu
                        temp1Str = append(temp1Str,temp1(i-1)) ; 
                        temp2Str = append(temp2Str,temp2(j-1)) ;  

                        %kontrol satırları
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

            %eşleşme matriste tersten okuma yaptığımız için okunan
            %değerlere flip işlemi
          
            temp1Str = flip(cell2mat(temp1Str)) ;
            temp2Str = flip(cell2mat(temp2Str)) ;
            fprintf('temp1Str : %s \n', temp1Str) ; 
            fprintf('temp2Str : %s \n', temp2Str) ; 
           
            % temp1 
            %başalngıç eşleşen indis değerleri üzerinden eşleşmeyen
            %kısımların yazdırılması
            if (i ~= 2)
                
                unMatchedStr = temp1(1:i-1);
                %eşleşen ve eşleşmeyen kısımların append olayı
                unMatchedStr = append(unMatchedStr, temp1Str) ; 
                
                if (row ~= NL + 1)
                    
                    unMatchedStr = append(unMatchedStr, temp1(row:length(temp1)));
        
                else 
                end
            
                
            elseif ( i == 2 )
                %eşleşmenin dizilerin  1. indisinden başlama durumu
                unMatchedStr = temp1(1) ;

                unMatchedStr = append(unMatchedStr, temp1Str) ;
                
                unMatchedStr = append(unMatchedStr, temp1(row:length(temp1)));
               
            end 
            


            %temp2str
            if (j ~= 2)
                
                unMatchedStr2 = temp2(1:j-1);
                unMatchedStr2 = append(unMatchedStr2, temp2Str) ; 
                
                if (column ~= NL + 1)
                    
                    unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
        
                
                

                else 
                end
                

            elseif ( j == 2 )
                unMatchedStr2 = temp2(1) ;

                unMatchedStr2 = append(unMatchedStr2, temp2Str) ;
                
                unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
                
            end  

            %dosyaya hizalanmış matrisleri yazdırma
            if (i > j)
                
                fprintf(file,'Hizalanmış 1 : %s\n', unMatchedStr);
                fprintf(file,'Hizalanmış 2 : ');
                %layout oluştuma boşluk atama 
                fprintf(file,repmat(' ',1,i-j));
                fprintf(file,'%s\n',unMatchedStr2);
                fprintf(file,'------------------------------------------------\n');
            
            elseif(j>i)
                
                fprintf(file,'Hizalanmış 1 : ');
                %layout oluştuma boşluk atama 
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



% READ COMPLEMENT 
fprintf(file, 'SEQUENCE COMPLEMENT KARŞILAŞTIRMASI\n');
fprintf('SEQUENCE COMPLEMENT KARŞILAŞTIRMASI\n');

for x = 1 : length(complementArray)
    temp1 = cell2mat(complementArray(x)) ;
    for y = x+1 : length(complementArray)
        temp2 = cell2mat(sequenceArray(y)) ;
        fprintf("Temp1 : %s \n", temp1) ; 
        fprintf("Temp2 : %s \n", temp2) ;

        
        readComplementMatrix = zeros(NL + 1, NL + 1) ; 
        %skor belirleme
       
        for k = 2 : length(temp1) + 1
            for d = 2 : length(temp2) + 1
                if (temp1(k - 1) == temp2( d - 1 ))
                    score2 = match ; 
                else
                    score2 = mismatch * -1 ; 
                end
                readComplementMatrix(k,d) = max([readComplementMatrix(k-1, d) - indel, readComplementMatrix(k-1, d-1) + score2, readComplementMatrix(k, d-1) - indel ]);

            end % for d 
            %maksimum belirleme
            [bestCompOverlap,column] = max(max(readComplementMatrix)) ;
            [~,row] = max(readComplementMatrix(:,column));
            
        end % for k

      


        %score yazdırma

        if bestCompOverlap >= score 
            overlapScoreMatrix(y,x) = bestCompOverlap ; 
            
    
            
            temp1Str = "" ;
            temp2Str = "" ;
        
    
            i = row ; 
            j = column ; 
            while ( i ~= 1  &&  j~=1 )
                 
                    %eşleşme durumlarını matriste tersten giderek yazdırma
                    
                    if readComplementMatrix(i, j-1) == readComplementMatrix(i,j) + indel
                        % SOL
                        temp1Str = append(temp1Str, '-') ; 
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                      
                        if j ~= 1
                            j = j-1 ; 
                        end
                   
                    elseif readComplementMatrix(i-1, j) == readComplementMatrix(i,j) + indel
                        % AŞAĞI YUKARI
                        temp1Str = append(temp1Str, temp1(i-1)) ; 
                        temp2Str = append(temp2Str, '-') ;  
                       
                        
                         if i ~= 1
                            i = i-1 ;
                            
                         end
                        
                    
                    elseif readComplementMatrix(i-1, j-1) == readComplementMatrix(i,j) - match
                        % ÇAPRAZ
                        temp1Str = append(temp1Str,temp1(i-1)) ;
                        temp2Str = append(temp2Str,temp2(j-1)) ;  
                        
                        
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
          %eşleşenlerin ters çevrilme durumu
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
    
              
            elseif ( i == 2 )
    
                unMatchedStr = temp1(1) ;
    
                unMatchedStr = append(unMatchedStr, temp1Str) ;
                
                unMatchedStr = append(unMatchedStr, temp1(row:length(temp1)));
                
            end %% temp1 str kontrol if i 
            
    
    
            %temp2str
            if (j ~= 2)
                
                unMatchedStr2 = temp2(1:j-1);
                unMatchedStr2 = append(unMatchedStr2, temp2Str) ; 
                
                if (column ~= NL + 1)
                    
                    unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
        
                
                
    
                else 
                end
               
    
            elseif ( j == 2 )
                unMatchedStr2 = temp2(1) ;
    
                unMatchedStr2 = append(unMatchedStr2, temp2Str) ;
                
                unMatchedStr2 = append(unMatchedStr2, temp2(column:length(temp2)));
               
            end %% temp1 str kontrol if i 
    
            %dosyaya hizalı yazdırma olayı
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


%dosyaya matris yazdırma işlemi
writematrix(overlapScoreMatrix,'bestoverlapscore.txt','Delimiter','tab');
     
fclose(file);

















