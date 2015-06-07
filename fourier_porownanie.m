%% porównanie sta³ego rozk³adu z jitterem dla transformaty fouriera
w_zadane.f_start=10;%wartoœci zadane
w_zadane.f_stop=400000;
w_zadane.krok=10;
w_zadane.l_probek=64;
karta_pomiarowa.f=4e7;
karta_pomiarowa.T=1/karta_pomiarowa.f;
w_pomocnicze.rozmiar_M=w_zadane.f_stop/w_zadane.krok;%przygotowanie macierzy zeby przyspieszyc czas
%równomierny rozk³ad próbek
const_rozklad.dft.fft=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
const_rozklad.dft.modul=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
const_rozklad.dft.faza=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
const_rozklad.x=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
const_rozklad.y=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
%równomierny rozk³ad próbek
const_rozklad_2.dft.fft=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek+1);
const_rozklad_2.dft.modul=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek+1);
const_rozklad_2.dft.faza=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek+1);
const_rozklad_2.x=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek+1);
const_rozklad_2.y=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek+1);
%nierównomierny rozk³ad próbek
jitter.dft.fft=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
jitter.dft.modul=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
jitter.dft.faza=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
jitter.x=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
jitter.y=zeros(w_pomocnicze.rozmiar_M,w_zadane.l_probek);
tic
j=0;
for i=w_zadane.f_start:w_zadane.krok:w_zadane.f_stop
    j=j+1;
    const_rozklad.f_zad(j)=i;
    const_rozklad.T_zad(j)=1/i; 
    const_rozklad.preskaler(j)=floor(const_rozklad.T_zad(j)/(w_zadane.l_probek*karta_pomiarowa.T));
    const_rozklad.f_pomiaru_probki(j)=karta_pomiarowa.f/const_rozklad.preskaler(j);
    const_rozklad.T_pomiaru_probki(j)=1/const_rozklad.f_pomiaru_probki(j);
    const_rozklad.x(j,:)=0:const_rozklad.T_pomiaru_probki(j):(w_zadane.l_probek-1)*const_rozklad.T_pomiaru_probki(j);
    const_rozklad.y(j,:)=sin(2*pi*i*const_rozklad.x(j,:));%+sin(8*pi*i*const_rozklad.x(j,:));
    const_rozklad.dft.fft(j,:)=fft(const_rozklad.y(j,:));
    const_rozklad.dft.modul(j,:)=abs(const_rozklad.dft.fft(j,:)/w_zadane.l_probek*2);
    const_rozklad.dft.faza(j,:)=angle(const_rozklad.dft.fft(j,:));
    const_rozklad.blad_T(j)=(const_rozklad.T_zad(j)-w_zadane.l_probek*const_rozklad.T_pomiaru_probki(j))/const_rozklad.T_zad(j)*100;
    %
    jitter.f_zad(j)=i; 
    jitter.T_zad(j)=const_rozklad.T_zad(j);
    jitter.l_T(j)=floor(jitter.T_zad(j)/karta_pomiarowa.T);
    jitter.x(j,:)=floor(linspace(0,jitter.l_T(j),w_zadane.l_probek))*karta_pomiarowa.T;
    jitter.y(j,:)=sin(2*pi*i*jitter.x(j,:));%+sin(8*pi*i*jitter.x(j,:));
    jitter.dft.fft(j,:)=fft(jitter.y(j,:));
    jitter.dft.modul(j,:)=abs(jitter.dft.fft(j,:)/w_zadane.l_probek*2);
    jitter.dft.faza(j,:)=angle(jitter.dft.fft(j,:));
    jitter.blad_T(j)=(jitter.T_zad(j)-jitter.x(j,end))/jitter.T_zad(j)*100;
    %
    const_rozklad_2.f_zad(j)=i;
    const_rozklad_2.T_zad(j)=1/i; 
    const_rozklad_2.preskaler(j)=floor(const_rozklad_2.T_zad(j)/((w_zadane.l_probek-1)*karta_pomiarowa.T));
    const_rozklad_2.f_pomiaru_probki(j)=karta_pomiarowa.f/const_rozklad_2.preskaler(j);
    const_rozklad_2.T_pomiaru_probki(j)=1/const_rozklad_2.f_pomiaru_probki(j);
    const_rozklad_2.x(j,:)=0:const_rozklad_2.T_pomiaru_probki(j):(w_zadane.l_probek)*const_rozklad_2.T_pomiaru_probki(j);
    const_rozklad_2.y(j,:)=sin(2*pi*i*const_rozklad_2.x(j,:));
    const_rozklad_2.dft.fft(j,:)=fft(const_rozklad_2.y(j,:));
    const_rozklad_2.dft.modul(j,:)=abs(const_rozklad_2.dft.fft(j,:)/(w_zadane.l_probek+1)*2);
    const_rozklad_2.dft.faza(j,:)=angle(const_rozklad_2.dft.fft(j,:));
    const_rozklad_2.blad_T(j)=(const_rozklad_2.T_zad(j)-(w_zadane.l_probek+1)*const_rozklad_2.T_pomiaru_probki(j))/const_rozklad_2.T_zad(j)*100;
end
toc