#PRZYGOTOWANIE BIBILOTEK
library(fastcluster)
library(igraph)
library(genie)
library(RSpectra)
library(speccalt)


#UPLOAD DANYCH
PATH <- file.path(getwd(), "pd2-zbiory-benchmarkowe")
datalist <- list.files(PATH, "data\\.gz$", recursive=TRUE)
labelslist <- list.files(PATH, "labels0\\.gz$", recursive=TRUE)
#listy plikow ze wszystkim danymi i etykietami

path_fcps <- file.path(getwd(), "pd2-zbiory-benchmarkowe", "fcps")
path_graves <- file.path(getwd(), "pd2-zbiory-benchmarkowe", "graves")
path_other <- file.path(getwd(), "pd2-zbiory-benchmarkowe", "other")
path_sipu <- file.path(getwd(), "pd2-zbiory-benchmarkowe", "sipu")
path_wut <- file.path(getwd(), "pd2-zbiory-benchmarkowe", "wut")


list_fcps <- list.files(path_fcps, "data\\.gz$", recursive=TRUE)
list_graves <- list.files(path_graves, "data\\.gz$", recursive=TRUE)
list_other <- list.files(path_other, "data\\.gz$", recursive=TRUE)
list_sipu <- list.files(path_sipu, "data\\.gz$", recursive=TRUE)
list_wut <- list.files(path_wut, "data\\.gz$", recursive=TRUE)


labels_fcps <- list.files(path_fcps, "labels0\\.gz$", recursive=TRUE)
labels_graves <- list.files(path_graves, "labels0\\.gz$", recursive=TRUE)
labels_other <- list.files(path_other, "labels0\\.gz$", recursive=TRUE)
labels_sipu <- list.files(path_sipu, "labels0\\.gz$", recursive=TRUE)
labels_wut <- list.files(path_wut, "labels0\\.gz$", recursive=TRUE)
#zaladowanie list plikow z danymi i etykietami podzielonymi na foldery

#FUNKCJE GENERUJACE
hclust_analyze <- function(data_set, labels_set, N, filename){
  # N - ilosc powtorzen wywolan funkcji
  labels_number <- max(labels_set)
  #ustalanie ilosci podzialow na podstawie danych etykiet
  for(i in c("average", "single", "complete", "mcquitty", "ward.D", "ward.D2", "centroid", "median")){
    for(j in 1:N){
      HC <- hclust(dist(data_set), method = i)
      #generowanie podzialu dla kazdej z metod
      div <- cutree(HC, labels_number)
      #zapisywanie podzialu do zmiennej
      output <- matrix(c("hclust", i, as.numeric(dendextend::FM_index(div, labels_set)), mclust::adjustedRandIndex(div, labels_set)), 
                       1, byrow = TRUE)
      write.table(output, file = filename, append = TRUE, col.names = FALSE, row.names = FALSE)
      #wypisanie do pliku wygenerowanej porcji danych
    }
  }
  
}


genie_analyze <- function(data_set, labels_set, N, filename){ 
  #zmienne odpowiednio 1. punkty 2. etykietu 3. ilosc powtorzen 4. nazwa pliku do zapisu
  labels_number <- max(labels_set)
  #ustalanie ilosci podzialow na podstawie danych etykiet
  for(i in c(0.1, 0.3, 0.5, 0.7, 0.9)){
    #wykonanie algorytmu dla 5 wsoplczynnikow
    for(j in 1:N){
      HC <- hclust2(dist(data_set), thresholdGini = i)
      div <- cutree(HC, labels_number)
      output <- matrix(c("Genie", paste("thresholdGini", i), as.numeric(dendextend::FM_index(div, labels_set)), mclust::adjustedRandIndex(div, labels_set)),
                       1, byrow = TRUE)
      write.table(output, file = filename, append = TRUE, col.names = FALSE, row.names = FALSE)
      #wypisanie do pliku wygenerowanej porcji danych
    }
  }
  
}


spectral_clustering_analyze <- function(data_set, labels_set, N, filename, M){
  #M - zmienna do algorytmu
  labels_number <- max(labels_set)
  #ustalanie ilosci podzialow na podstawie danych etykiet
  for(i in 1:M){
    div <- spectral_clustering_0(data_set, M, labels_number)
    #funkcja generujaca "podstawe" podzialu
    for(j in 1:N){
      div <- kmeans(div, labels_number)$cluster
      #element zmienny w generowaniu podzialu, koncowy etap generowania podzialu
      output <- matrix(c("Own Algorythm", paste("M", i), as.numeric(dendextend::FM_index(div, labels_set)), mclust::adjustedRandIndex(div, labels_set)),
                   1, byrow = TRUE)
      write.table(output, file = filename, append = TRUE, col.names = FALSE, row.names = FALSE)
      #wypisanie do pliku wygenerowanej porcji danych
    }
  }
}

speccalt_analyze <- function(data_set, labels_set, N, filename){
  labels_number <- max(labels_set)
  #ustalanie ilosci podzialow na podstawie danych etykiet
  for(i in 1:N){
    kernel_matrix <- local.rbfdot(data_set)
    #generowanie jadra macierzy (funkcja z biblioteki speccalt)
    div <- speccalt(kernel_matrix, labels_number)
    #generowanie podzialu
    output <- matrix(c("speccalt", "-", as.numeric(dendextend::FM_index(div, labels_set)), mclust::adjustedRandIndex(div, labels_set)),
                     1, byrow = TRUE)
    write.table(output, file = filename, append = TRUE, col.names = FALSE, row.names = FALSE)
    #wypisanie do pliku wygenerowanej porcji danych
  }
}

#PARAMETRYZOWANIE SPECTRAL CLUSTERING I POWTARZANIA TESTOW

M <- 10 #zmienna do spectral_clustering
N <- 1000#ilosc powtorzen wywolania kazdej funkcji

#GENEROWANIE CALOSCIOWE
generating_all <- function(){
  for(i in 1:length(datalist)){
    data_set <- read.table(file.path(PATH, datalist[i]))
    labels_set <- read.table(file.path(PATH, labelslist[i]))$V1
  
    hclust_analyze(data_set, labels_set, N, paste("s", i))
    genie_analyze(data_set, labels_set, N, paste("s", i)) 
    spectral_clustering_analyze(data_set, labels_set, N, paste("s", i), M)
    speccalt_analyze(data_set, labels_set, N, "spec")
    #wlasciwie nie ma dobrego powodu, zeby robic to inaczej niz pozostale funkcje
    #z jakiegos w danej chwili tak bylo mi wygodniej zrzucic to do jednego pliku
  }
}



#GENEROWANIE CZESCIOWE 
generating_part <- function(){
  for(i in 3:length(list_fcps)){
    data_set <- read.table(file.path(path_fcps, list_fcps[i]))
    labels_set <- read.table(file.path(path_fcps, labels_fcps[i]))$V1
  
    hclust_analyze(data_set, labels_set, N, paste("fcps", i))
    genie_analyze(data_set, labels_set, N, paste("fcps", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("fcps", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("fcps", i))
  
  }

  for(i in 1:length(list_graves)){
    data_set <- read.table(file.path(path_graves, list_graves[i]))
    labels_set <- read.table(file.path(path_graves, labels_graves[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("graves", i))
    genie_analyze(data_set, labels_set, N, paste("graves", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("graves", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("graves", i))
  }

  for(i in 1:length(list_other)){
    data_set <- read.table(file.path(path_other, list_other[i]))
    labels_set <- read.table(file.path(path_other, labels_other[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("other", i))
    genie_analyze(data_set, labels_set, N, paste("other", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("other", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("other", i))
  }

  for(i in 1:length(list_sipu)){
    data_set <- read.table(file.path(path_sipu, list_sipu[i]))
    labels_set <- read.table(file.path(path_sipu, labels_sipu[i]))$V1
  
    hclust_analyze(data_set, labels_set, N, paste("sipu", i))
    genie_analyze(data_set, labels_set, N, paste("sipu", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("sipu", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("sipu", i))
  }

  for(i in 1:length(list_wut)){
    data_set <- read.table(file.path(path_wut, list_wut[i]))
    labels_set <- read.table(file.path(path_wut, labels_wut[i]))$V1
  
    hclust_analyze(data_set, labels_set, N, paste("wut", i))
    genie_analyze(data_set, labels_set, N, paste("wut", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("wut", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("wut", i))
  }


}

#GENEROWANIE CALOSCIOWE PO STANDARYZACJI
generating_all_norm <- function(){
  for(i in 1:length(datalist)){
    data_set <- scale(read.table(file.path(PATH, datalist[i])))
    labels_set <- read.table(file.path(PATH, labelslist[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("sn", i))
    genie_analyze(data_set, labels_set, N, paste("sn", i)) 
    spectral_clustering_analyze(data_set, labels_set, N, paste("sn", i), M)
  }
  speccalt_analyze(data_set, labels_set, N, "specn")
}



#GENEROWANIE CZESCIOWE PO STANDARYZACJI
generating_part_norm <- function(){
  for(i in 1:length(list_fcps)){
    data_set <- scale(read.table(file.path(path_fcps, list_fcps[i])))
    labels_set <- read.table(file.path(path_fcps, labels_fcps[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("fcpsn", i))
    genie_analyze(data_set, labels_set, N, paste("fcpsn", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("fcpsn", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("fcpsn", i))
    
  }
  
  for(i in 1:length(list_graves)){
    data_set <- scale(read.table(file.path(path_graves, list_graves[i])))
    labels_set <- read.table(file.path(path_graves, labels_graves[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("gravesn", i))
    genie_analyze(data_set, labels_set, N, paste("gravesn", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("gravesn", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("gravesn", i))
  }
  
  for(i in 1:length(list_other)){
    data_set <- scale(read.table(file.path(path_other, list_other[i])))
    labels_set <- read.table(file.path(path_other, labels_other[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("othern", i))
    genie_analyze(data_set, labels_set, N, paste("othern", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("othern", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("othern", i))
  }
  
  for(i in 1:length(list_sipu)){
    data_set <- scale(read.table(file.path(path_sipu, list_sipu[i])))
    labels_set <- read.table(file.path(path_sipu, labels_sipu[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("sipun", i))
    genie_analyze(data_set, labels_set, N, paste("sipun", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("sipun", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("sipun", i))
  }
  
  for(i in 1:length(list_wut)){
    data_set <- scale(read.table(file.path(path_wut, list_wut[i])))
    labels_set <- read.table(file.path(path_wut, labels_wut[i]))$V1
    
    hclust_analyze(data_set, labels_set, N, paste("wutn", i))
    genie_analyze(data_set, labels_set, N, paste("wutn", i))
    spectral_clustering_analyze(data_set, labels_set, N, paste("wutn", i), M)
    speccalt_analyze(data_set, labels_set, N, paste("wutn", i))
  }
  
  
}








