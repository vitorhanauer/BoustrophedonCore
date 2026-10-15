# Visão Geral

O intuito do projeto é fazer uma implementação do algoritmo Boustrophedon para percorrer um poligono e entregar uma array de pontos. Essa implementação deverá ser feita em C. Essa será a lógica core, porém, essa implementação será utilizadas por outras linguagens de programação, ou seja, as outras linguagens farão um wrapper dessa implementação. Essa implementação também precisa conseguir se chamada por outros códigos em C.

O processos do algortimo de dividem em alguns segmentos:
1. Dividir a área recebida em polígonos menores continuos.
2. Tratar os polígonos menores como nós de um grafo e encontrar o menor caminho entre esses nós.
3. Gerar o caminho de pontos em cada nó
4. Unir os caminhos de cada nó o, o final de um nó se liga ao inicio do outro nó.

Cada uma dessas etapas é será descrita posteriormente

# Assinatura da função

A implementação deverá ter uma função de ponto de entrada, ou seja, uma única função que outros código irão chamar. Isso não significa que todo o código deve ficar em uma única função. Abaixo o pseudo código da assinatura da função.

```
generateWaypoints(Ponto Inicial, Ponto Final, Vertices da Area, Vertices dos Bloqueios, Largura de atuação, Sobreposição);
```
**Ponto Inicial**: Representado por dois valores, [10,10] ou {x: 10, y: 10}. Será de onde o agente autonomo irá começar o percurso.
**Ponto Final**: Representado por dois valores, [10,10] ou {x: 10, y: 10}. Será de onde o agente autonomo irá terminar o percurso.
**Vertices da Area**: Representado por vários valores, [[10,10],[10,9],...] ou [{x: 10, y:10}, {x:10, y: 9}]. Representa os vertices do poligono a ser percorrido.
**Vertices dos Bloqueios**: Representado por vários valores, [[[10,10],[10,9],...], [[10,10],[10,9],...]] ou [[{x: 10, y:10}, {x:10, y: 9}], [{x: 10, y:10}, {x:10, y: 9}]]. Representa os vertices do obstáculos.
**Largura de atuação**: Representado por um único valor, 10. Representa a largura total da visão do agente autonomo, valor da direita + valor da esquerda.
**Sobreposição**: Representado por um único valor, 0 a 100. Representa o quanto as larguras de visão devem se sobrepor.

OBS: A distância entre as linhas é calculada através da seguinte formular: distância = largura(1 - sobreposição).

# Práticas de codificação

- Utilizar convetional commits e sempre commitar quando tiver terminado a menor unidade de uma atividade funcional.
- Utilizar nomes de variáveis, funções etc de forma que fiquem legíveis.
- Utilizar TDD para codificar primeiro um teste funcional e depois criar a funcionalidade.

# Limitação

Para o cenário de teste será utilizado uma configuração de 512MB de RAM e 1Ghz de clock. 
Os polígonos de testes devem ser de:
- 10x10 utilizando larguras de atuação de 1,5.
- 100x100 utilizando larguras de atuação de 1,5,10,50.
- 1_000x1_000 utilizando larguras de atuação de 1,5,10,50.
- 10_000x10_000 utilizando larguras de atuação de 1,5,10,50.
