#!/bin/python3

'''
    This program provides a way of implementing DES algorithm in python3 
    without use of framework/library outside Python standard libraries.
    The program provides all functionality for each of  the 16 Feistel
    rounds, generation of the 16 48-bits subkeys and finally
    crypting/decrypting plain/crypted text respectively. 
    Copyright (C) 2016  Man-gel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

class des:
    'Clase que implementa la funcionalidad del algoritmo DES'
    __name__= 'des'
    __file__='des.py'
    
    __TAM_BLQ=64
    __N_RONDAS=16
    
    __PC1_TABLA = { 0:56,  1:48,  2:40,  3:32,  4:24,  5:16,  6:8,   7:0,
                    8:57,  9:49, 10:41, 11:33, 12:25, 13:17, 14:9,  15:1,
                   16:58, 17:50, 18:42, 19:34, 20:26, 21:18, 22:10, 23:2,
                   24:59, 25:51, 26:43, 27:35, 28:62, 29:54, 30:46, 31:38,
                   32:30, 33:22, 34:14, 35:6,  36:61, 37:53, 38:45, 39:37,
                   40:29, 41:21, 42:13, 43:5,  44:60, 45:52, 46:44, 47:36,
                   48:28, 49:20, 50:12, 51:4,  52:27, 53:19, 54:11, 55:3  }
    
    __PC2_TABLA = { 0:13,  1:16,  2:10,  3:23,  4:0,   5:4,   6:2,   7:27,
                    8:14,  9:5,  10:20, 11:9,  12:22, 13:18, 14:11, 15:3,
                   16:25, 17:7,  18:15, 19:6,  20:26, 21:19, 22:12, 23:1,
                   24:40, 25:51, 26:30, 27:36, 28:46, 29:54, 30:29, 31:39,
                   32:50, 33:44, 34:32, 35:47, 36:43, 37:48, 38:38, 39:55,
                   40:33, 41:52, 42:45, 43:41, 44:49, 45:35, 46:28, 47:31  }
    
    __IP_TABLA={ 0:57,  1:49,  2:41,  3:33,  4:25,  5:17,  6:9,   7:1,
                 8:59,  9:51, 10:43, 11:35, 12:27, 13:19, 14:11, 15:3,
                16:61, 17:53, 18:45, 19:37, 20:29, 21:21, 22:13, 23:5,
                24:63, 25:55, 26:47, 27:39, 28:31, 29:23, 30:15, 31:7,
                32:56, 33:48, 34:40, 35:32, 36:24, 37:16, 38:8,  39:0,
                40:58, 41:50, 42:42, 43:34, 44:26, 45:18, 46:10, 47:2,
                48:60, 49:52, 50:44, 51:36, 52:28, 53:20, 54:12, 55:4,
                56:62, 57:54, 58:46, 59:38, 60:30, 61:22, 62:14, 63:6 }
    
    __IP_INV_TABLA={ 0:39,  1:7,   2:47,  3:15,  4:55,  5:23,  6:63,  7:31,
                     8:38,  9:6,  10:46, 11:14, 12:54, 13:22, 14:62, 15:30,
                    16:37, 17:5,  18:45, 19:13, 20:53, 21:21, 22:61, 23:29,
                    24:36, 25:4,  26:44, 27:12, 28:52, 29:20, 30:60, 31:28,
                    32:35, 33:3,  34:43, 35:11, 36:51, 37:19, 38:59, 39:27,
                    40:34, 41:2,  42:42, 43:10, 44:50, 45:18, 46:58, 47:26,
                    48:33, 49:1,  50:41, 51:9,  52:49, 53:17, 54:57, 55:25,
                    56:32, 57:0,  58:40, 59:8,  60:48, 61:16, 62:56, 63:24 }
    
    __E_TABLA={3:4, 7:8, 11:12, 15:16, 19:20, 23:24, 27:28}
    
    __S_BOX={ 0:{ 0:(14,4,13,1,2,15,11,8,3,10,6,12,5,9,0,7),   #S_BOX 1
                  1:(0,15,7,4,14,2,13,1,10,6,12,11,9,5,3,8),
                  2:(4,1,14,8,13,6,2,11,15,12,9,7,3,10,5,0),
                  3:(15,12,8,2,4,9,1,7,5,11,3,14,10,0,6,13) },     
    
              6:{ 0:(15,1,8,14,6,11,3,4,9,7,2,13,12,0,5,10),   #S_BOX 2
                  1:(3,13,4,7,15,2,8,14,12,0,1,10,6,9,11,5),
                  2:(0,14,7,11,10,4,13,1,5,8,12,6,9,3,2,15),
                  3:(13,8,10,1,3,15,4,2,11,6,7,12,0,5,14,9) },
    
             12:{ 0:(10,0,9,14,6,3,15,5,1,13,12,7,11,4,2,8),   #S_BOX 3
                  1:(13,7,0,9,3,4,6,10,2,8,5,14,12,11,15,1),
                  2:(13,6,4,9,8,15,3,0,11,1,2,12,5,10,14,7),
                  3:(1,10,13,0,6,9,8,7,4,15,14,3,11,5,2,12) },
      
             18:{ 0:(7,13,14,3,0,6,9,10,1,2,8,5,11,12,4,15),   #S_BOX 4
                  1:(13,8,11,5,6,15,0,3,4,7,2,12,1,10,14,9),
                  2:(10,6,9,0,12,11,7,13,15,1,3,14,5,2,8,4),
                  3:(3,15,0,6,10,1,13,8,9,4,5,11,12,7,2,14) },
      
             24:{ 0:(2,12,4,1,7,10,11,6,8,5,3,15,13,0,14,9),   #S_BOX 5
                  1:(14,11,2,12,4,7,13,1,5,0,15,10,3,9,8,6),
                  2:(4,2,1,11,10,13,7,8,15,9,12,5,6,3,0,14),
                  3:(11,8,12,7,1,14,2,13,6,15,0,9,10,4,5,3) },
      
             30:{ 0:(12,1,10,15,9,2,6,8,0,13,3,4,14,7,5,11),   #S_BOX 6
                  1:(10,15,4,2,7,12,9,5,6,1,13,14,0,11,3,8),
                  2:(9,14,15,5,2,8,12,3,7,0,4,10,1,13,11,6),
                  3:(4,3,2,12,9,5,15,10,11,14,1,7,6,0,8,13) },
      
             36:{ 0:(4,11,2,14,15,0,8,13,3,12,9,7,5,10,6,1),   #S_BOX 7
                  1:(13,0,11,7,4,9,1,10,14,3,5,12,2,15,8,6),
                  2:(1,4,11,13,12,3,7,14,10,15,6,8,0,5,9,2),
                  3:(6,11,13,8,1,4,10,7,9,5,0,15,14,2,3,12)  },
      
             42:{ 0:(13,2,8,4,6,15,11,1,10,9,3,14,5,0,12,7),   #S_BOX 8
                  1:(1,15,13,8,10,3,7,4,12,5,6,11,0,14,9,2),
                  2:(7,11,4,1,9,12,14,2,0,6,10,13,15,3,5,8),
                  3:(2,1,14,7,4,10,8,13,15,12,9,0,3,5,6,11)  }}

             
    __PERM_TABLA = { 0:15,  1:6,   2:19,  3:20,  4:28,  5:11,  6:27,  7:16,
                     8:0,   9:14, 10:22, 11:25, 12:4,  13:17, 14:30, 15:9,
                    16:1,  17:7,  18:23, 19:13, 20:31, 21:26, 22:2,  23:8,
                    24:18, 25:12, 26:29, 27:5,  28:21, 29:10, 30:3,  31:24  }
        
    #s ->(string) mensaje + bits de relleno. Regresa list de strings, c/u = c/letra del msg en su representación binaria.
    def __stringTobits(self, s): 
        res = []
        mask='00000000'
        for c in s:
            bUnmasked= bin( ord(c) )[2:]
            bMasked=mask[ len(bUnmasked): ] + bUnmasked
            res.append(bMasked)
        return res

    '''
    Según PKCS#5 si faltan n bytes se debe comletar con el valor n, n veces el octeto.
    Las posibilidades van desde 1 a 7 que en decimal van desde 49 a 55 (ASCII/UNICODE).
    Si n = 1 el bloque se completaría con '1' =>(49 en dec, 00110001 en bin), si n=2 el boque se completaría con '22' =>(5050 en dec, 0011001000110010 en bin), etc.
    '''
    
    #msg ->(string) mensaje del usuario. Regresa(entero) cant de caracteres que faltan para que len(msg) sea múltiplo de 8    
    def __faltaRelleno(self, msg):
        lets = 0
        letsXblq = self.__TAM_BLQ / 8
        for l in msg:
            lets += 1
        if lets % letsXblq == 0: #Regresa 0 si es múltiplo del tamaño de bloque y no necesita rellenarse, SINO regresa x (1<x<7) 
            return 0
        else:
            rell = letsXblq - (lets % letsXblq)
            return rell

    '''
    msg -> mensaje del usuario(string). faltantes -> cant faltante(entero). Regresa msg(como string) concatenado a chars de relleno (ver PKCS#5)
    '''
    def __rellenar(self, msg,faltantes):
        if faltantes == 1:
            return msg+'1'
        elif faltantes == 2:
            return msg+'22'    
        elif faltantes == 3:
            return msg+'333'
        elif faltantes == 4:
            return msg+'4444'
        elif faltantes == 5:
            return msg+'55555'
        elif faltantes == 6:
            return msg+'666666'
        elif faltantes == 7:
            return msg+'7777777'

    #msg -> mensaje desencriptado(string).Regresa msg sin el patrón de relleno(como string) 
    def __desrellenar(self,msg):
        for i in range(8):
            if i != 0:
                if msg[len(msg)-1] == str(i) and msg[len(msg)-i] == str(i):                    
                    return msg[:(len(msg)-i)]
        return msg
    #part -> list de strings de len()=8 (c/u = 1 caracter del msg como octeto). Regresa list de strings con len()=64
    def __bloques64( self, part ): 
        bloque=''
        bits=0
        blq64 = []
        for oc in range( len(part) ):
            for bit in range( len(part[oc]) ):
                bloque += part[oc][bit]
                bits += 1
                if bits == self.__TAM_BLQ:
                    blq64.append(bloque)
                    bloque=''
                    bits = 0
        return tuple(blq64)

    def __calcularClaves(self, k): #k -> string o list de chars len()=64, key original. Regresa list con 16 keys (strings len()=48 c/u).
        kList= []
        kPerm = []
        izq=2
        for b in range( len(self.__PC1_TABLA) ):
            kPerm.append( k[ self.__PC1_TABLA[b] ] )
        C0 = kPerm[: int(len(kPerm)/2) ]
        D0 = kPerm[ int(len(kPerm)/2) :]
        for r in range(self.__N_RONDAS):
            if r == 0 or r == 1 or r == 8 or r == 15:
                izq = 1
            else:
                izq = 2
            rotC = C0[:izq]
            rotD = D0[:izq]
            for i in range(28-izq):
                C0[i] = C0[i+izq]
                D0[i] = D0[i+izq]
            del C0[28-izq:]
            del D0[28-izq:]
            C0 = C0 + rotC
            D0 = D0 + rotD
            kConc = C0
            kConc.extend(D0)
            Ki=[]
            for i in range( len(self.__PC2_TABLA) ):
                Ki.append( kConc[self.__PC2_TABLA[i]] )
            kList.append(Ki)
        return tuple(kList)

    def __pInicial(self, blq): #blq -> string de len()=64 con '0'|'1'. Regresa list de enteros con 0 | 1 de len()=64
        blqPerm = []
        for b in range( len(blq) ):
            blqPerm.append(int( blq[ self.__IP_TABLA[b] ] ))        
        return tuple(blqPerm)

    def __pFinal(self, blq): #blq -> string de len()=64 con '0'|'1'. Regresa list de enteros con 0 | 1 de len()=64
        blqPerm = []
        for b in range( len(blq) ):
            blqPerm.append(int( blq[ self.__IP_INV_TABLA[b] ] ))        
        return tuple(blqPerm)

    def __expansion(self, blq): #blq ->  list = 1 bloque de 32 bits(como enteros), Regresa list de 48 bits (como enteros)    
        blqExp=[]
        blqExp.append( blq[-1] )
        for b in range( len(blq) ):
            if b in self.__E_TABLA:
                blqExp.append( blq[b] )
                blqExp.append( blq[ self.__E_TABLA[b] ] )
                blqExp.append( blq[b] )
            else:
                blqExp.append( blq[b] )
        blqExp.append( blq[0] )
        return tuple(blqExp)

    #bits -> list de enteros(1|0). k -> list o string solo con chars '0'|'1'.Siempre con len(bits). Regresa list de enteros 0|1.
    def __xor(self, bits, k): 
        _xor = []
        for i in range(len(bits)):
            if bits[i] == 0 and str(k[i]) == '1' or bits[i] == 1 and str(k[i]) == '0':
                _xor.append(1)
            else:
                _xor.append(0)
        return tuple(_xor)

    def __sustitucion(self, bits): #bits -> list de enteros 0 | 1 con len()=48. Regresa list con len()=32 con enteros 0| 1
        bitSust = []
        mask = '0000'
        cont = 0
        while cont < 48:
            nbit = ''
            col=''
            row=''
            r=0
            c=0
            blq = bits[cont:cont+6]
            row=( str(blq[0])+str(blq[-1]) )
            if row != '00':
                r = int(row,2)
            col=(str(blq[1]) + str(blq[2]) + str(blq[3]) + str(blq[4]) )
            if col != '0000':
                c = int(col,2)
            nbit = str( bin(self.__S_BOX[cont][r][c]) )
            nbit = nbit[2:]
            nbit = mask[len(nbit):] + nbit        
            for nb in nbit:
                bitSust.append( int(nb) )        
            cont += 6        
        return tuple(bitSust)

    def __permutacion(self, blq): #blq -> list de len()=32 con enteros 1 | 0. Regresa una estructura igual permutada.
        bPerm = []
        for b in range( len(blq) ):
            bPerm.append( blq[ self.__PERM_TABLA[b] ] )
        return tuple(bPerm)

    #blqL, blqR -> list len()=32 c/u (enteros).Regresa string len()=64 con bloques intercambiados
    def __intercambiarBloques(self, blqL, blqR): 
        blqExch = ''
        for r in blqR:
            blqExch += str(r)
        for l in blqL:
            blqExch += str(l)    
        return blqExch


    ''' 
    m -> bloque de mensaje, string len()=64. k -> claves de c/ronda (list de 16 strings len()=48 - Ki -). Regresa bloque encriptado (list de enteros 0|1 len()=64 ).
    '''
    def __encriptar(self, m, k):    
        ronda = 0
        li=[]                                 #li es para mantener Li-1
        ri=[]                                 #ri es para mantener Ri-1
        Ri=[]                                 #Ri se usará en las iteraciones siguientes a la primera
        Li=[]                                 #Li se usará en las iteraciones siguientes a la primera
        pi = self.__pInicial(m)               #permutar con IP y obtener blq 64bits
        L0 = pi[:int(len(pi)/2)]              #dividir en blq inicial L0 de 32bits
        R0 = pi[int(len(pi)/2):]              #dividir en blq inicial R0 de 32bits
        while ronda < self.__N_RONDAS:        #i desde 1 a 16
            if ronda > 0:
                li = Li                       #guardar en li -> Li-1 (Li ahora es Li-1)                    ecuación 7.4
                ri = Ri                       #guardar en ri -> Ri-1 (Ri ahora es Ri-1)                    ecuación 7.4
            else:
                li = L0                       #guardar en li -> Li-1 (L0 al inicio es Li-1)                ecuación 7.4
                ri = R0                       #guardar en ri -> Ri-1 (R0 al inicio es Ri-1)                ecuación 7.4
                del L0
                del R0
            e = self.__expansion(ri)          #expansión de R0|Ri de 32bits a 48bits                       2.4.1
            x = self.__xor( e, k[ronda] )     #aplicar xor entre e y k de 48bits                           2.4.2
            s = self.__sustitucion(x)         #sustituir blq de 48bits por 32bits usando sboxes1-8         2.4.3 y 2.4.4
            f = self.__permutacion(s)         #permutar blq 32bits y obtener f(Ri-1, Ki)                   2.4.5
            Ri = self.__xor( li, f )          #guardar Ri -> Li-1 xor f(Ri-1, Ki)(xor 32bits)(li es Li-1)  ecuación 7.5
            Li = ri                           #aplicar Li = Ri-1 (el intercambio para la siguiente ronda)  Fig. 7.9
            ronda += 1                        #pasar a la siguiente ronda y seguir con otro key
        i = self.__intercambiarBloques(Li,Ri) #intercambiar últimos bloques y unirlos en uno de 64bits
        return self.__pFinal(i)               #permutar con IP inversa y obtener encriptación de 64bits

    #hx ->(string) msg|key en hexadecimal. Regresa list de strings, c/u = c/par letras de hx(1octeto) en representación binaria.
    def __hexTobin(self, hx):
        toHex = []
        mask='00000000'
        x=0
        while x < len(hx):
            h=''
            B = str( bin( int( (str(hx[x:x+2])),16 ) ) )[2:]
            h = mask[ len(B):]+B
            x += 2
            toHex.append(h)
        return toHex

    def __binToStr(self, b):#b ->  list de enteros len()=64. Regresa string equivalente a las letras de cada octeto de la lista.
        w=''
        oc=''
        ch = ''
        for i in range( len(b) ):        
            oc += str( b[i] )
            if i > 0 and (i+1)%8 == 0:
                ch=str( int(oc,2) )
                w += chr( int(ch) )
                oc=''                
        return w

    def __keyHextoBin(self, k): #k -> string con la clave en hexadecimal. Regresa string len()=64 (representación binaria de la clave).
        kBin = ''
        mask = '00000000'
        i=0
        while i < len(k):
            unmsk = str( bin( int(k[i:i+2],16) ) )[2:]
            kBin += ( mask[len(unmsk):]  + unmsk )
            i += 2
        return kBin
    
    ''' 
    m -> bloque de mensaje encriptado, string len()=64. k -> claves de c/ronda (list de 16 strings len()=48 - Ki -). Regresa bloque desencriptado (list de enteros 0|1 len()=64 ).
    '''    
    def __desencriptar(self, m,k):
        ronda = 15
        li=[]                                 #li es para mantener Li-1
        ri=[]                                 #ri es para mantener Ri-1
        Ri=[]                                 #Ri se usará en las iteraciones siguientes a la primera
        Li=[]                                 #Li se usará en las iteraciones siguientes a la primera
        pi = self.__pInicial(m)               #permutar con IP y obtener blq 64bits
        L0 = pi[:int(len(pi)/2)]              #dividir en blq inicial L0 de 32bits
        R0 = pi[int(len(pi)/2):]              #dividir en blq inicial R0 de 32bits
        while ronda > -1:                     #i desde 16 a 1
            if ronda < self.__N_RONDAS -1 :
                li = Li                       #guardar en li -> Li-1 (Li ahora es Li-1)                    ecuación 7.4
                ri = Ri                       #guardar en ri -> Ri-1 (Ri ahora es Ri-1)                    ecuación 7.4
            else:
                li = L0                       #guardar en li -> Li-1 (L0 al inicio es Li-1)                ecuación 7.4
                ri = R0                       #guardar en ri -> Ri-1 (R0 al inicio es Ri-1)                ecuación 7.4
                del L0
                del R0
            e = self.__expansion(ri)          #expansión de R0|Ri 32bits  a 48bits                         2.4.1
            x = self.__xor( e, k[ronda] )     #aplicar xor a blq 48bits con k 48bits                       2.4.2
            s = self.__sustitucion(x)         #sustituir blq 48bits por 32bits con respectiva sbox         2.4.3 y 2.4.4
            f = self.__permutacion(s)         #permutar blq 32bits a otro igual                            2.4.5
            Ri = self.__xor( li, f )          #guardar Ri -> Li-1 xor f(Ri-1, Ki). (li es Li-1)            ecuación 7.5
            Li = ri                           #aplicar Li = Ri-1 (el intercambio para la siguiente ronda)  Fig. 7.9
            ronda -= 1                        #pasar a la siguiente ronda y seguir con otro key
        i = self.__intercambiarBloques(Li,Ri) #intercambiar últimos bloques y unirlos en uno de 64bits
        return self.__pFinal(i)               #permutar con IP inversa y obtener encriptación de 64bits



    '''MÉTODOS CON TODA LA FUNCIONALIDAD PARA IMPLEMENTAR EL ALGORITMO, ES LO QUE NECESITAS PARA ENCRIPTAR CON DES'''
    ''' TODOS LOS METODOS ANTERIORES SON PARA USO INTERNO DE LA CLASE'''
    '''
    word -> string = password. Regresa -1 si word no es valida o tupla. Donde: 
    tupla[0]= word/password <- hexadecimal (Es para poder compartirla y que otro pueda desencriptar, solo opcionalmente)
    tupla[1]=list con 16 keys (strings len()=48 c/u).
    '''
    def genKeys(self,word):
        if len(word) == 16:
            k = self.__keyHextoBin(word)
            return ( word,self.__calcularClaves(k))         
        elif len(word) == 8:
            kBin = ""
            kHex = ""
            wBin = self.__stringTobits(word)
            for s in wBin:
                kBin += s
                kHex += str( hex( int(s,2) ) )[2:]
            return ( kHex,self.__calcularClaves(kBin))            
        else:
            return -1
        
    #m -> mensaje a encriptar(string). keys-> list con 16 keys (strings len()=48 c/u). Regresa mensaje encriptado en hexadecimal (string).
    def crypt(self,m,keys):
        if(len(keys) != 16):
            return "Error: 16 keys invalidas"
        mEnc = []
        h = ''
        f = self.__faltaRelleno(m)
        if f != 0:
            m=self.__rellenar(m,f)
        mbin=self.__stringTobits(m)
        blq64 = self.__bloques64(mbin)
        for blq in blq64:
            h += hex( int(blq,2) )[2:]
            mEnc.append( self.__encriptar(blq,keys) )
        cr=''
        mcry=''
        for bq in mEnc:
            for b in bq:
                cr += str(b)
            mcry += str( hex( int(cr,2) ) )[2:]
            cr=''
        return mcry
    
    #m -> mensaje a desencriptar(string). keys-> list con 16 keys (strings len()=48 c/u). Regresa mensaje desencriptado (string-texto plano).
    def decrypt(self,m,keys):
        if(len(keys) != 16):
            return "Error: 16 keys invalidas"
        mDes = []
        mbin = []
        mbin=self.__hexTobin(m)
        blq64 = self.__bloques64(mbin)
        for blq in blq64:
            mDes.append( self.__desencriptar(blq,keys) )
        msg=''
        for blqB in mDes:
            msg += self.__binToStr(blqB)
        return self.__desrellenar(msg)
