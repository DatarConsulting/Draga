
p4_arithmetic/code/arithmetic_filter.o:     file format elf64-x86-64


Disassembly of section .text:

0000000000000000 <weight>:
       0:	53                   	push   %rbx
       1:	48 89 f9             	mov    %rdi,%rcx
       4:	85 f6                	test   %esi,%esi
       6:	0f 84 8c 01 00 00    	je     198 <weight+0x198>
       c:	48 b8 00 00 00 00 00 	movabs $0x3000000000000,%rax
      13:	00 03 00 
      16:	4c 8d 8c 3f 00 00 00 	lea    0x5000000(%rdi,%rdi,1),%r9
      1d:	05 
      1e:	41 ba 00 00 00 03    	mov    $0x3000000,%r10d
      24:	48 99                	cqto
      26:	49 f7 f9             	idiv   %r9
      29:	49 89 c0             	mov    %rax,%r8
      2c:	49 8d 81 00 00 00 fa 	lea    -0x6000000(%r9),%rax
      33:	49 8d b1 00 00 00 06 	lea    0x6000000(%r9),%rsi
      3a:	48 c1 e0 18          	shl    $0x18,%rax
      3e:	48 99                	cqto
      40:	48 f7 fe             	idiv   %rsi
      43:	49 89 c3             	mov    %rax,%r11
      46:	48 0f af c0          	imul   %rax,%rax
      4a:	4c 89 db             	mov    %r11,%rbx
      4d:	48 c1 f8 18          	sar    $0x18,%rax
      51:	48 0f af d8          	imul   %rax,%rbx
      55:	48 85 db             	test   %rbx,%rbx
      58:	48 8d b3 ff ff ff 00 	lea    0xffffff(%rbx),%rsi
      5f:	48 0f 49 f3          	cmovns %rbx,%rsi
      63:	48 c1 fe 18          	sar    $0x18,%rsi
      67:	48 0f af f0          	imul   %rax,%rsi
      6b:	48 b8 67 66 66 66 66 	movabs $0x6666666666666667,%rax
      72:	66 66 66 
      75:	48 f7 ee             	imul   %rsi
      78:	48 c1 fe 3f          	sar    $0x3f,%rsi
      7c:	48 b8 ab aa aa aa aa 	movabs $0x2aaaaaaaaaaaaaab,%rax
      83:	aa aa 2a 
      86:	48 c1 fa 19          	sar    $0x19,%rdx
      8a:	48 89 d7             	mov    %rdx,%rdi
      8d:	48 f7 eb             	imul   %rbx
      90:	48 c1 fb 3f          	sar    $0x3f,%rbx
      94:	48 b8 09 e1 d1 c6 11 	movabs $0x29f16b11c6d1e109,%rax
      9b:	6b f1 29 
      9e:	48 29 f7             	sub    %rsi,%rdi
      a1:	48 c1 fa 17          	sar    $0x17,%rdx
      a5:	48 29 da             	sub    %rbx,%rdx
      a8:	5b                   	pop    %rbx
      a9:	4c 01 da             	add    %r11,%rdx
      ac:	48 01 d7             	add    %rdx,%rdi
      af:	48 8d b4 3f 17 72 b1 	lea    0xb17217(%rdi,%rdi,1),%rsi
      b6:	00 
      b7:	48 69 f6 d6 07 00 00 	imul   $0x7d6,%rsi,%rsi
      be:	48 f7 ee             	imul   %rsi
      c1:	48 8d 81 00 00 00 01 	lea    0x1000000(%rcx),%rax
      c8:	48 c1 fe 3f          	sar    $0x3f,%rsi
      cc:	49 0f af c2          	imul   %r10,%rax
      d0:	48 c1 fa 0f          	sar    $0xf,%rdx
      d4:	48 89 d7             	mov    %rdx,%rdi
      d7:	48 99                	cqto
      d9:	49 f7 f9             	idiv   %r9
      dc:	49 b9 43 78 b4 71 c4 	movabs $0xa7c5ac471b47843,%r9
      e3:	5a 7c 0a 
      e6:	48 29 f7             	sub    %rsi,%rdi
      e9:	48 8d b0 00 00 80 ff 	lea    -0x800000(%rax),%rsi
      f0:	48 b8 56 55 55 55 55 	movabs $0x5555555555555556,%rax
      f7:	55 55 55 
      fa:	48 f7 ee             	imul   %rsi
      fd:	48 c1 fe 3f          	sar    $0x3f,%rsi
     101:	48 89 d1             	mov    %rdx,%rcx
     104:	48 29 f1             	sub    %rsi,%rcx
     107:	48 0f af c9          	imul   %rcx,%rcx
     10b:	48 c1 f9 18          	sar    $0x18,%rcx
     10f:	48 89 ca             	mov    %rcx,%rdx
     112:	48 0f af d1          	imul   %rcx,%rdx
     116:	48 69 c9 f9 84 01 00 	imul   $0x184f9,%rcx,%rcx
     11d:	48 c1 fa 18          	sar    $0x18,%rdx
     121:	48 69 d2 fb ea 01 00 	imul   $0x1eafb,%rdx,%rdx
     128:	48 c1 ea 05          	shr    $0x5,%rdx
     12c:	48 89 d0             	mov    %rdx,%rax
     12f:	49 f7 e1             	mul    %r9
     132:	48 89 d6             	mov    %rdx,%rsi
     135:	48 c1 ee 07          	shr    $0x7,%rsi
     139:	48 c1 e9 05          	shr    $0x5,%rcx
     13d:	48 89 c8             	mov    %rcx,%rax
     140:	49 f7 e1             	mul    %r9
     143:	48 89 f8             	mov    %rdi,%rax
     146:	48 0f af c7          	imul   %rdi,%rax
     14a:	48 c1 ea 07          	shr    $0x7,%rdx
     14e:	48 c1 f8 19          	sar    $0x19,%rax
     152:	48 8d 94 16 00 00 00 	lea    0x1000000(%rsi,%rdx,1),%rdx
     159:	01 
     15a:	48 29 f8             	sub    %rdi,%rax
     15d:	48 8d 88 00 00 00 01 	lea    0x1000000(%rax),%rcx
     164:	49 0f af c8          	imul   %r8,%rcx
     168:	48 85 c9             	test   %rcx,%rcx
     16b:	48 8d 81 ff ff ff 00 	lea    0xffffff(%rcx),%rax
     172:	48 0f 49 c1          	cmovns %rcx,%rax
     176:	48 c1 f8 18          	sar    $0x18,%rax
     17a:	48 0f af d0          	imul   %rax,%rdx
     17e:	48 85 d2             	test   %rdx,%rdx
     181:	48 8d 82 ff ff ff 00 	lea    0xffffff(%rdx),%rax
     188:	48 0f 49 c2          	cmovns %rdx,%rax
     18c:	48 c1 f8 18          	sar    $0x18,%rax
     190:	c3                   	ret
     191:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
     198:	41 b9 00 00 80 03    	mov    $0x3800000,%r9d
     19e:	41 ba 00 00 00 01    	mov    $0x1000000,%r10d
     1a4:	48 b8 00 00 00 00 00 	movabs $0x3000000000000,%rax
     1ab:	00 03 00 
     1ae:	49 29 f9             	sub    %rdi,%r9
     1b1:	48 99                	cqto
     1b3:	4d 01 c9             	add    %r9,%r9
     1b6:	49 f7 f9             	idiv   %r9
     1b9:	49 89 c0             	mov    %rax,%r8
     1bc:	e9 6b fe ff ff       	jmp    2c <weight+0x2c>
     1c1:	66 66 2e 0f 1f 84 00 	data16 cs nopw 0x0(%rax,%rax,1)
     1c8:	00 00 00 00 
     1cc:	0f 1f 40 00          	nopl   0x0(%rax)

00000000000001d0 <output.isra.0>:
     1d0:	41 57                	push   %r15
     1d2:	49 89 d7             	mov    %rdx,%r15
     1d5:	48 89 f8             	mov    %rdi,%rax
     1d8:	31 d2                	xor    %edx,%edx
     1da:	41 56                	push   %r14
     1dc:	4d 89 fe             	mov    %r15,%r14
     1df:	41 55                	push   %r13
     1e1:	45 31 ed             	xor    %r13d,%r13d
     1e4:	41 54                	push   %r12
     1e6:	55                   	push   %rbp
     1e7:	4c 89 cd             	mov    %r9,%rbp
     1ea:	49 89 f1             	mov    %rsi,%r9
     1ed:	53                   	push   %rbx
     1ee:	48 89 cb             	mov    %rcx,%rbx
     1f1:	48 89 f9             	mov    %rdi,%rcx
     1f4:	4c 89 7c 24 c0       	mov    %r15,-0x40(%rsp)
     1f9:	45 31 ff             	xor    %r15d,%r15d
     1fc:	4c 01 f0             	add    %r14,%rax
     1ff:	49 89 c2             	mov    %rax,%r10
     202:	48 89 54 24 d0       	mov    %rdx,-0x30(%rsp)
     207:	48 89 c8             	mov    %rcx,%rax
     20a:	4c 11 fa             	adc    %r15,%rdx
     20d:	49 89 d3             	mov    %rdx,%r11
     210:	48 89 74 24 a8       	mov    %rsi,-0x58(%rsp)
     215:	48 89 de             	mov    %rbx,%rsi
     218:	48 03 74 24 a8       	add    -0x58(%rsp),%rsi
     21d:	48 89 7c 24 c8       	mov    %rdi,-0x38(%rsp)
     222:	4c 89 ef             	mov    %r13,%rdi
     225:	48 c7 44 24 b0 00 00 	movq   $0x0,-0x50(%rsp)
     22c:	00 00 
     22e:	48 13 7c 24 b0       	adc    -0x50(%rsp),%rdi
     233:	48 f7 e3             	mul    %rbx
     236:	4c 89 6c 24 e0       	mov    %r13,-0x20(%rsp)
     23b:	48 89 4c 24 e8       	mov    %rcx,-0x18(%rsp)
     240:	44 89 c1             	mov    %r8d,%ecx
     243:	49 89 c4             	mov    %rax,%r12
     246:	49 89 d5             	mov    %rdx,%r13
     249:	4c 89 c8             	mov    %r9,%rax
     24c:	4c 89 4c 24 f0       	mov    %r9,-0x10(%rsp)
     251:	48 f7 64 24 c0       	mulq   -0x40(%rsp)
     256:	48 89 5c 24 d8       	mov    %rbx,-0x28(%rsp)
     25b:	49 29 c4             	sub    %rax,%r12
     25e:	b8 01 00 00 00       	mov    $0x1,%eax
     263:	49 19 d5             	sbb    %rdx,%r13
     266:	48 d3 e0             	shl    %cl,%rax
     269:	4c 89 ea             	mov    %r13,%rdx
     26c:	48 0f af d0          	imul   %rax,%rdx
     270:	48 89 d1             	mov    %rdx,%rcx
     273:	49 f7 e4             	mul    %r12
     276:	48 01 ca             	add    %rcx,%rdx
     279:	49 89 c0             	mov    %rax,%r8
     27c:	49 89 d1             	mov    %rdx,%r9
     27f:	49 c1 e0 05          	shl    $0x5,%r8
     283:	49 0f a4 c1 05       	shld   $0x5,%rax,%r9
     288:	49 01 c0             	add    %rax,%r8
     28b:	4c 89 d0             	mov    %r10,%rax
     28e:	49 11 d1             	adc    %rdx,%r9
     291:	4c 89 da             	mov    %r11,%rdx
     294:	4c 01 d0             	add    %r10,%rax
     297:	4d 0f a4 c1 04       	shld   $0x4,%r8,%r9
     29c:	49 c1 e0 04          	shl    $0x4,%r8
     2a0:	4c 0f a4 d2 01       	shld   $0x1,%r10,%rdx
     2a5:	4d 89 c4             	mov    %r8,%r12
     2a8:	4c 01 d0             	add    %r10,%rax
     2ab:	49 89 f0             	mov    %rsi,%r8
     2ae:	4c 11 da             	adc    %r11,%rdx
     2b1:	4d 89 cd             	mov    %r9,%r13
     2b4:	49 c1 e0 02          	shl    $0x2,%r8
     2b8:	49 89 f9             	mov    %rdi,%r9
     2bb:	49 0f a4 f1 02       	shld   $0x2,%rsi,%r9
     2c0:	49 01 f0             	add    %rsi,%r8
     2c3:	49 11 f9             	adc    %rdi,%r9
     2c6:	4c 01 c0             	add    %r8,%rax
     2c9:	4d 89 d0             	mov    %r10,%r8
     2cc:	4c 11 ca             	adc    %r9,%rdx
     2cf:	4d 89 d9             	mov    %r11,%r9
     2d2:	49 c1 e0 02          	shl    $0x2,%r8
     2d6:	4d 0f a4 d1 02       	shld   $0x2,%r10,%r9
     2db:	4d 01 d0             	add    %r10,%r8
     2de:	49 89 f2             	mov    %rsi,%r10
     2e1:	4d 11 d9             	adc    %r11,%r9
     2e4:	49 89 fb             	mov    %rdi,%r11
     2e7:	49 01 f2             	add    %rsi,%r10
     2ea:	49 0f a4 f3 01       	shld   $0x1,%rsi,%r11
     2ef:	49 01 f2             	add    %rsi,%r10
     2f2:	48 89 c6             	mov    %rax,%rsi
     2f5:	48 89 d0             	mov    %rdx,%rax
     2f8:	49 11 fb             	adc    %rdi,%r11
     2fb:	4d 01 d0             	add    %r10,%r8
     2fe:	4d 11 d9             	adc    %r11,%r9
     301:	49 0f af c0          	imul   %r8,%rax
     305:	4c 89 ca             	mov    %r9,%rdx
     308:	48 0f af d6          	imul   %rsi,%rdx
     30c:	48 89 c1             	mov    %rax,%rcx
     30f:	48 89 f0             	mov    %rsi,%rax
     312:	48 01 d1             	add    %rdx,%rcx
     315:	49 f7 e0             	mul    %r8
     318:	48 01 ca             	add    %rcx,%rdx
     31b:	48 89 c6             	mov    %rax,%rsi
     31e:	48 01 c6             	add    %rax,%rsi
     321:	48 89 d7             	mov    %rdx,%rdi
     324:	48 0f a4 c7 01       	shld   $0x1,%rax,%rdi
     329:	48 01 c6             	add    %rax,%rsi
     32c:	48 11 d7             	adc    %rdx,%rdi
     32f:	48 0f a4 f7 04       	shld   $0x4,%rsi,%rdi
     334:	48 c1 e6 04          	shl    $0x4,%rsi
     338:	48 01 c6             	add    %rax,%rsi
     33b:	48 11 d7             	adc    %rdx,%rdi
     33e:	4c 39 e6             	cmp    %r12,%rsi
     341:	48 89 f8             	mov    %rdi,%rax
     344:	4c 19 e8             	sbb    %r13,%rax
     347:	0f 82 33 01 00 00    	jb     480 <output.isra.0+0x2b0>
     34d:	48 8b 7c 24 e8       	mov    -0x18(%rsp),%rdi
     352:	b9 03 00 00 00       	mov    $0x3,%ecx
     357:	4c 8b 4c 24 f0       	mov    -0x10(%rsp),%r9
     35c:	be 05 00 00 00       	mov    $0x5,%esi
     361:	48 89 f8             	mov    %rdi,%rax
     364:	48 f7 e1             	mul    %rcx
     367:	49 89 c4             	mov    %rax,%r12
     36a:	4c 89 c8             	mov    %r9,%rax
     36d:	49 89 d5             	mov    %rdx,%r13
     370:	48 f7 e6             	mul    %rsi
     373:	49 01 c4             	add    %rax,%r12
     376:	48 89 f8             	mov    %rdi,%rax
     379:	48 8b 7c 24 c0       	mov    -0x40(%rsp),%rdi
     37e:	49 11 d5             	adc    %rdx,%r13
     381:	48 f7 e6             	mul    %rsi
     384:	49 89 c6             	mov    %rax,%r14
     387:	4c 89 c8             	mov    %r9,%rax
     38a:	49 89 d7             	mov    %rdx,%r15
     38d:	48 f7 e1             	mul    %rcx
     390:	49 01 c6             	add    %rax,%r14
     393:	48 89 f8             	mov    %rdi,%rax
     396:	49 11 d7             	adc    %rdx,%r15
     399:	48 f7 e1             	mul    %rcx
     39c:	49 89 c0             	mov    %rax,%r8
     39f:	48 89 d8             	mov    %rbx,%rax
     3a2:	49 89 d1             	mov    %rdx,%r9
     3a5:	48 f7 e6             	mul    %rsi
     3a8:	49 01 c0             	add    %rax,%r8
     3ab:	48 89 f8             	mov    %rdi,%rax
     3ae:	49 11 d1             	adc    %rdx,%r9
     3b1:	4d 01 e0             	add    %r12,%r8
     3b4:	4d 11 e9             	adc    %r13,%r9
     3b7:	48 f7 e6             	mul    %rsi
     3ba:	48 89 c6             	mov    %rax,%rsi
     3bd:	48 89 d8             	mov    %rbx,%rax
     3c0:	48 89 d7             	mov    %rdx,%rdi
     3c3:	48 f7 e1             	mul    %rcx
     3c6:	4c 89 c9             	mov    %r9,%rcx
     3c9:	48 01 c6             	add    %rax,%rsi
     3cc:	48 11 d7             	adc    %rdx,%rdi
     3cf:	4c 01 f6             	add    %r14,%rsi
     3d2:	4c 11 ff             	adc    %r15,%rdi
     3d5:	48 0f af ce          	imul   %rsi,%rcx
     3d9:	48 89 74 24 a8       	mov    %rsi,-0x58(%rsp)
     3de:	48 89 f8             	mov    %rdi,%rax
     3e1:	48 89 7c 24 b0       	mov    %rdi,-0x50(%rsp)
     3e6:	49 0f af c0          	imul   %r8,%rax
     3ea:	48 01 c1             	add    %rax,%rcx
     3ed:	4c 89 c0             	mov    %r8,%rax
     3f0:	48 f7 e6             	mul    %rsi
     3f3:	48 01 ca             	add    %rcx,%rdx
     3f6:	48 89 c1             	mov    %rax,%rcx
     3f9:	49 89 c2             	mov    %rax,%r10
     3fc:	48 89 d3             	mov    %rdx,%rbx
     3ff:	49 89 d3             	mov    %rdx,%r11
     402:	49 c1 e2 02          	shl    $0x2,%r10
     406:	48 0f a4 cb 01       	shld   $0x1,%rcx,%rbx
     40b:	48 01 c9             	add    %rcx,%rcx
     40e:	49 0f a4 c3 02       	shld   $0x2,%rax,%r11
     413:	48 01 c1             	add    %rax,%rcx
     416:	48 89 f8             	mov    %rdi,%rax
     419:	48 11 d3             	adc    %rdx,%rbx
     41c:	49 0f af c4          	imul   %r12,%rax
     420:	4c 89 ea             	mov    %r13,%rdx
     423:	48 0f af d6          	imul   %rsi,%rdx
     427:	4c 89 fe             	mov    %r15,%rsi
     42a:	48 8d 3c 02          	lea    (%rdx,%rax,1),%rdi
     42e:	4c 89 e0             	mov    %r12,%rax
     431:	48 f7 64 24 a8       	mulq   -0x58(%rsp)
     436:	48 01 fa             	add    %rdi,%rdx
     439:	48 29 c1             	sub    %rax,%rcx
     43c:	4c 89 c8             	mov    %r9,%rax
     43f:	48 19 d3             	sbb    %rdx,%rbx
     442:	49 0f af c6          	imul   %r14,%rax
     446:	49 0f af f0          	imul   %r8,%rsi
     44a:	48 01 c6             	add    %rax,%rsi
     44d:	4c 89 f0             	mov    %r14,%rax
     450:	49 f7 e0             	mul    %r8
     453:	48 01 f2             	add    %rsi,%rdx
     456:	48 29 c1             	sub    %rax,%rcx
     459:	48 19 d3             	sbb    %rdx,%rbx
     45c:	48 89 4d 00          	mov    %rcx,0x0(%rbp)
     460:	48 89 5d 08          	mov    %rbx,0x8(%rbp)
     464:	4c 89 55 10          	mov    %r10,0x10(%rbp)
     468:	4c 89 5d 18          	mov    %r11,0x18(%rbp)
     46c:	5b                   	pop    %rbx
     46d:	5d                   	pop    %rbp
     46e:	41 5c                	pop    %r12
     470:	41 5d                	pop    %r13
     472:	41 5e                	pop    %r14
     474:	41 5f                	pop    %r15
     476:	c3                   	ret
     477:	66 0f 1f 84 00 00 00 	nopw   0x0(%rax,%rax,1)
     47e:	00 00 
     480:	48 8b 74 24 d8       	mov    -0x28(%rsp),%rsi
     485:	48 8b 7c 24 e0       	mov    -0x20(%rsp),%rdi
     48a:	48 8b 44 24 c8       	mov    -0x38(%rsp),%rax
     48f:	48 8b 54 24 d0       	mov    -0x30(%rsp),%rdx
     494:	48 03 44 24 a8       	add    -0x58(%rsp),%rax
     499:	48 13 54 24 b0       	adc    -0x50(%rsp),%rdx
     49e:	4c 01 f6             	add    %r14,%rsi
     4a1:	4c 11 ff             	adc    %r15,%rdi
     4a4:	48 01 c6             	add    %rax,%rsi
     4a7:	48 11 d7             	adc    %rdx,%rdi
     4aa:	48 89 f1             	mov    %rsi,%rcx
     4ad:	49 89 f2             	mov    %rsi,%r10
     4b0:	49 89 fb             	mov    %rdi,%r11
     4b3:	48 89 fb             	mov    %rdi,%rbx
     4b6:	48 01 f1             	add    %rsi,%rcx
     4b9:	49 c1 e2 02          	shl    $0x2,%r10
     4bd:	48 0f a4 f3 01       	shld   $0x1,%rsi,%rbx
     4c2:	49 0f a4 f3 02       	shld   $0x2,%rsi,%r11
     4c7:	48 01 f1             	add    %rsi,%rcx
     4ca:	48 11 fb             	adc    %rdi,%rbx
     4cd:	48 0f a4 c2 01       	shld   $0x1,%rax,%rdx
     4d2:	48 01 c0             	add    %rax,%rax
     4d5:	48 29 c1             	sub    %rax,%rcx
     4d8:	48 19 d3             	sbb    %rdx,%rbx
     4db:	e9 7c ff ff ff       	jmp    45c <output.isra.0+0x28c>

00000000000004e0 <ac_branch>:
     4e0:	f3 0f 1e fa          	endbr64
     4e4:	55                   	push   %rbp
     4e5:	48 01 fe             	add    %rdi,%rsi
     4e8:	b8 18 00 00 00       	mov    $0x18,%eax
     4ed:	53                   	push   %rbx
     4ee:	31 db                	xor    %ebx,%ebx
     4f0:	48 01 ff             	add    %rdi,%rdi
     4f3:	48 01 db             	add    %rbx,%rbx
     4f6:	48 39 f7             	cmp    %rsi,%rdi
     4f9:	72 07                	jb     502 <ac_branch+0x22>
     4fb:	48 29 f7             	sub    %rsi,%rdi
     4fe:	48 83 c3 01          	add    $0x1,%rbx
     502:	83 e8 01             	sub    $0x1,%eax
     505:	75 e9                	jne    4f0 <ac_branch+0x10>
     507:	48 89 df             	mov    %rbx,%rdi
     50a:	31 f6                	xor    %esi,%esi
     50c:	e8 ef fa ff ff       	call   0 <weight>
     511:	48 89 df             	mov    %rbx,%rdi
     514:	be 01 00 00 00       	mov    $0x1,%esi
     519:	48 89 c5             	mov    %rax,%rbp
     51c:	e8 df fa ff ff       	call   0 <weight>
     521:	5b                   	pop    %rbx
     522:	48 89 c2             	mov    %rax,%rdx
     525:	48 89 e8             	mov    %rbp,%rax
     528:	48 c1 e0 10          	shl    $0x10,%rax
     52c:	48 01 d5             	add    %rdx,%rbp
     52f:	31 d2                	xor    %edx,%edx
     531:	48 f7 f5             	div    %rbp
     534:	5d                   	pop    %rbp
     535:	c3                   	ret
     536:	66 2e 0f 1f 84 00 00 	cs nopw 0x0(%rax,%rax,1)
     53d:	00 00 00 

0000000000000540 <decode>:
     540:	41 57                	push   %r15
     542:	41 56                	push   %r14
     544:	41 55                	push   %r13
     546:	41 54                	push   %r12
     548:	55                   	push   %rbp
     549:	53                   	push   %rbx
     54a:	48 89 d3             	mov    %rdx,%rbx
     54d:	8d 57 fc             	lea    -0x4(%rdi),%edx
     550:	48 83 ec 58          	sub    $0x58,%rsp
     554:	48 89 74 24 28       	mov    %rsi,0x28(%rsp)
     559:	83 fa 1a             	cmp    $0x1a,%edx
     55c:	0f 87 40 03 00 00    	ja     8a2 <decode+0x362>
     562:	69 cf 8b 8a 01 00    	imul   $0x18a8b,%edi,%ecx
     568:	ba 01 00 00 00       	mov    $0x1,%edx
     56d:	48 89 d0             	mov    %rdx,%rax
     570:	81 c1 3f 0d 03 00    	add    $0x30d3f,%ecx
     576:	c1 e9 06             	shr    $0x6,%ecx
     579:	48 69 c9 63 2d 3e 05 	imul   $0x53e2d63,%rcx,%rcx
     580:	48 c1 e9 26          	shr    $0x26,%rcx
     584:	83 c1 03             	add    $0x3,%ecx
     587:	48 d3 e0             	shl    %cl,%rax
     58a:	48 39 c6             	cmp    %rax,%rsi
     58d:	0f 83 0f 03 00 00    	jae    8a2 <decode+0x362>
     593:	48 85 f6             	test   %rsi,%rsi
     596:	0f 84 e0 02 00 00    	je     87c <decode+0x33c>
     59c:	89 f9                	mov    %edi,%ecx
     59e:	c7 44 24 44 00 00 00 	movl   $0x0,0x44(%rsp)
     5a5:	00 
     5a6:	31 ff                	xor    %edi,%edi
     5a8:	45 31 ff             	xor    %r15d,%r15d
     5ab:	48 d3 e2             	shl    %cl,%rdx
     5ae:	48 89 5c 24 48       	mov    %rbx,0x48(%rsp)
     5b3:	48 89 c1             	mov    %rax,%rcx
     5b6:	bd 01 00 00 00       	mov    $0x1,%ebp
     5bb:	48 89 54 24 30       	mov    %rdx,0x30(%rsp)
     5c0:	41 be 01 00 00 00    	mov    $0x1,%r14d
     5c6:	48 89 fb             	mov    %rdi,%rbx
     5c9:	48 c7 44 24 38 00 00 	movq   $0x0,0x38(%rsp)
     5d0:	00 00 
     5d2:	48 c7 44 24 20 01 00 	movq   $0x1,0x20(%rsp)
     5d9:	00 00 
     5db:	48 c7 44 24 08 00 00 	movq   $0x0,0x8(%rsp)
     5e2:	00 00 
     5e4:	0f 1f 40 00          	nopl   0x0(%rax)
     5e8:	4d 89 fa             	mov    %r15,%r10
     5eb:	45 31 db             	xor    %r11d,%r11d
     5ee:	49 89 dc             	mov    %rbx,%r12
     5f1:	48 89 e8             	mov    %rbp,%rax
     5f4:	49 01 ea             	add    %rbp,%r10
     5f7:	49 83 d3 00          	adc    $0x0,%r11
     5fb:	45 31 ed             	xor    %r13d,%r13d
     5fe:	4d 01 f4             	add    %r14,%r12
     601:	49 83 d5 00          	adc    $0x0,%r13
     605:	49 f7 e6             	mul    %r14
     608:	49 89 c0             	mov    %rax,%r8
     60b:	48 89 d8             	mov    %rbx,%rax
     60e:	49 89 d1             	mov    %rdx,%r9
     611:	49 f7 e7             	mul    %r15
     614:	48 89 c6             	mov    %rax,%rsi
     617:	4c 89 c0             	mov    %r8,%rax
     61a:	4c 8b 44 24 30       	mov    0x30(%rsp),%r8
     61f:	48 89 d7             	mov    %rdx,%rdi
     622:	48 29 f0             	sub    %rsi,%rax
     625:	4c 89 ca             	mov    %r9,%rdx
     628:	4c 8b 4c 24 38       	mov    0x38(%rsp),%r9
     62d:	48 19 fa             	sbb    %rdi,%rdx
     630:	48 89 44 24 10       	mov    %rax,0x10(%rsp)
     635:	4c 89 c6             	mov    %r8,%rsi
     638:	48 8b 7c 24 10       	mov    0x10(%rsp),%rdi
     63d:	48 0f af f2          	imul   %rdx,%rsi
     641:	4c 89 ca             	mov    %r9,%rdx
     644:	48 0f af d7          	imul   %rdi,%rdx
     648:	48 89 f8             	mov    %rdi,%rax
     64b:	48 01 d6             	add    %rdx,%rsi
     64e:	49 f7 e0             	mul    %r8
     651:	48 01 f2             	add    %rsi,%rdx
     654:	48 89 c6             	mov    %rax,%rsi
     657:	48 89 d7             	mov    %rdx,%rdi
     65a:	48 c1 e6 05          	shl    $0x5,%rsi
     65e:	48 0f a4 c7 05       	shld   $0x5,%rax,%rdi
     663:	48 01 c6             	add    %rax,%rsi
     666:	4c 89 e0             	mov    %r12,%rax
     669:	48 11 d7             	adc    %rdx,%rdi
     66c:	4c 89 ea             	mov    %r13,%rdx
     66f:	48 0f a4 f7 04       	shld   $0x4,%rsi,%rdi
     674:	48 c1 e6 04          	shl    $0x4,%rsi
     678:	48 89 74 24 10       	mov    %rsi,0x10(%rsp)
     67d:	4c 89 d6             	mov    %r10,%rsi
     680:	48 89 7c 24 18       	mov    %rdi,0x18(%rsp)
     685:	4c 01 d6             	add    %r10,%rsi
     688:	4c 89 df             	mov    %r11,%rdi
     68b:	4c 0f a4 d7 01       	shld   $0x1,%r10,%rdi
     690:	4c 01 d6             	add    %r10,%rsi
     693:	4c 11 df             	adc    %r11,%rdi
     696:	48 c1 e0 02          	shl    $0x2,%rax
     69a:	49 89 f0             	mov    %rsi,%r8
     69d:	4c 89 e6             	mov    %r12,%rsi
     6a0:	4c 0f a4 e2 02       	shld   $0x2,%r12,%rdx
     6a5:	4c 01 e0             	add    %r12,%rax
     6a8:	49 89 f9             	mov    %rdi,%r9
     6ab:	4c 89 ef             	mov    %r13,%rdi
     6ae:	4c 11 ea             	adc    %r13,%rdx
     6b1:	49 01 c0             	add    %rax,%r8
     6b4:	4c 89 d0             	mov    %r10,%rax
     6b7:	49 11 d1             	adc    %rdx,%r9
     6ba:	4c 89 da             	mov    %r11,%rdx
     6bd:	48 c1 e0 02          	shl    $0x2,%rax
     6c1:	4c 0f a4 d2 02       	shld   $0x2,%r10,%rdx
     6c6:	4c 01 d0             	add    %r10,%rax
     6c9:	4c 11 da             	adc    %r11,%rdx
     6cc:	4c 01 e6             	add    %r12,%rsi
     6cf:	4c 0f a4 e7 01       	shld   $0x1,%r12,%rdi
     6d4:	4c 01 e6             	add    %r12,%rsi
     6d7:	4c 11 ef             	adc    %r13,%rdi
     6da:	48 01 f0             	add    %rsi,%rax
     6dd:	48 11 fa             	adc    %rdi,%rdx
     6e0:	4c 89 cf             	mov    %r9,%rdi
     6e3:	48 0f af f8          	imul   %rax,%rdi
     6e7:	49 0f af d0          	imul   %r8,%rdx
     6eb:	48 89 fe             	mov    %rdi,%rsi
     6ee:	48 01 d6             	add    %rdx,%rsi
     6f1:	49 f7 e0             	mul    %r8
     6f4:	48 01 f2             	add    %rsi,%rdx
     6f7:	48 89 c6             	mov    %rax,%rsi
     6fa:	48 89 d7             	mov    %rdx,%rdi
     6fd:	48 01 c6             	add    %rax,%rsi
     700:	48 0f a4 c7 01       	shld   $0x1,%rax,%rdi
     705:	48 01 c6             	add    %rax,%rsi
     708:	48 11 d7             	adc    %rdx,%rdi
     70b:	48 0f a4 f7 04       	shld   $0x4,%rsi,%rdi
     710:	48 c1 e6 04          	shl    $0x4,%rsi
     714:	48 01 c6             	add    %rax,%rsi
     717:	48 8b 44 24 10       	mov    0x10(%rsp),%rax
     71c:	48 11 d7             	adc    %rdx,%rdi
     71f:	48 8b 54 24 18       	mov    0x18(%rsp),%rdx
     724:	48 39 c6             	cmp    %rax,%rsi
     727:	48 19 d7             	sbb    %rdx,%rdi
     72a:	72 14                	jb     740 <decode+0x200>
     72c:	b8 fe ff ff ff       	mov    $0xfffffffe,%eax
     731:	48 83 c4 58          	add    $0x58,%rsp
     735:	5b                   	pop    %rbx
     736:	5d                   	pop    %rbp
     737:	41 5c                	pop    %r12
     739:	41 5d                	pop    %r13
     73b:	41 5e                	pop    %r14
     73d:	41 5f                	pop    %r15
     73f:	c3                   	ret
     740:	83 7c 24 44 14       	cmpl   $0x14,0x44(%rsp)
     745:	74 e5                	je     72c <decode+0x1ec>
     747:	49 8d 34 1e          	lea    (%r14,%rbx,1),%rsi
     74b:	49 8d 3c 2f          	lea    (%r15,%rbp,1),%rdi
     74f:	4c 8d 61 ff          	lea    -0x1(%rcx),%r12
     753:	e8 00 00 00 00       	call   758 <decode+0x218>
			754: R_X86_64_PLT32	ac_branch-0x4
     758:	89 c0                	mov    %eax,%eax
     75a:	49 0f af c4          	imul   %r12,%rax
     75e:	48 89 c1             	mov    %rax,%rcx
     761:	48 c1 e9 10          	shr    $0x10,%rcx
     765:	48 3d ff ff 00 00    	cmp    $0xffff,%rax
     76b:	0f 86 f7 00 00 00    	jbe    868 <decode+0x328>
     771:	4c 39 e1             	cmp    %r12,%rcx
     774:	0f 83 ee 00 00 00    	jae    868 <decode+0x328>
     77a:	48 8b 54 24 08       	mov    0x8(%rsp),%rdx
     77f:	48 83 c2 01          	add    $0x1,%rdx
     783:	48 8d 04 11          	lea    (%rcx,%rdx,1),%rax
     787:	48 39 44 24 28       	cmp    %rax,0x28(%rsp)
     78c:	48 89 44 24 08       	mov    %rax,0x8(%rsp)
     791:	0f 93 c0             	setae  %al
     794:	0f 83 96 00 00 00    	jae    830 <decode+0x2f0>
     79a:	48 8d 7c 2d 00       	lea    0x0(%rbp,%rbp,1),%rdi
     79f:	4c 8d 04 5b          	lea    (%rbx,%rbx,2),%r8
     7a3:	48 89 54 24 08       	mov    %rdx,0x8(%rsp)
     7a8:	48 8d 34 1f          	lea    (%rdi,%rbx,1),%rsi
     7ac:	48 01 ef             	add    %rbp,%rdi
     7af:	4f 8d 0c 76          	lea    (%r14,%r14,2),%r9
     7b3:	4a 8d 1c 47          	lea    (%rdi,%r8,2),%rbx
     7b7:	4f 8d 04 3f          	lea    (%r15,%r15,1),%r8
     7bb:	48 89 f5             	mov    %rsi,%rbp
     7be:	4b 8d 3c 30          	lea    (%r8,%r14,1),%rdi
     7c2:	4d 01 f8             	add    %r15,%r8
     7c5:	4f 8d 34 48          	lea    (%r8,%r9,2),%r14
     7c9:	49 89 ff             	mov    %rdi,%r15
     7cc:	48 8b 54 24 20       	mov    0x20(%rsp),%rdx
     7d1:	0f b6 c0             	movzbl %al,%eax
     7d4:	83 44 24 44 01       	addl   $0x1,0x44(%rsp)
     7d9:	48 8b 7c 24 08       	mov    0x8(%rsp),%rdi
     7de:	48 01 d2             	add    %rdx,%rdx
     7e1:	48 09 c2             	or     %rax,%rdx
     7e4:	48 89 54 24 20       	mov    %rdx,0x20(%rsp)
     7e9:	48 39 7c 24 28       	cmp    %rdi,0x28(%rsp)
     7ee:	0f 85 f4 fd ff ff    	jne    5e8 <decode+0xa8>
     7f4:	48 89 d8             	mov    %rbx,%rax
     7f7:	48 8b 5c 24 48       	mov    0x48(%rsp),%rbx
     7fc:	48 89 43 08          	mov    %rax,0x8(%rbx)
     800:	8b 44 24 44          	mov    0x44(%rsp),%eax
     804:	48 89 2b             	mov    %rbp,(%rbx)
     807:	89 43 20             	mov    %eax,0x20(%rbx)
     80a:	48 8b 44 24 20       	mov    0x20(%rsp),%rax
     80f:	4c 89 7b 10          	mov    %r15,0x10(%rbx)
     813:	48 89 43 28          	mov    %rax,0x28(%rbx)
     817:	31 c0                	xor    %eax,%eax
     819:	4c 89 73 18          	mov    %r14,0x18(%rbx)
     81d:	48 83 c4 58          	add    $0x58,%rsp
     821:	5b                   	pop    %rbx
     822:	5d                   	pop    %rbp
     823:	41 5c                	pop    %r12
     825:	41 5d                	pop    %r13
     827:	41 5e                	pop    %r14
     829:	41 5f                	pop    %r15
     82b:	c3                   	ret
     82c:	0f 1f 40 00          	nopl   0x0(%rax)
     830:	48 8d 14 1b          	lea    (%rbx,%rbx,1),%rdx
     834:	49 29 cc             	sub    %rcx,%r12
     837:	48 8d 7c 6d 00       	lea    0x0(%rbp,%rbp,2),%rdi
     83c:	48 8d 34 1a          	lea    (%rdx,%rbx,1),%rsi
     840:	48 8d 1c 2a          	lea    (%rdx,%rbp,1),%rbx
     844:	4c 89 e1             	mov    %r12,%rcx
     847:	4b 8d 14 36          	lea    (%r14,%r14,1),%rdx
     84b:	4f 8d 0c 7f          	lea    (%r15,%r15,2),%r9
     84f:	4e 8d 04 32          	lea    (%rdx,%r14,1),%r8
     853:	48 8d 2c 7e          	lea    (%rsi,%rdi,2),%rbp
     857:	4e 8d 34 3a          	lea    (%rdx,%r15,1),%r14
     85b:	4f 8d 3c 48          	lea    (%r8,%r9,2),%r15
     85f:	e9 68 ff ff ff       	jmp    7cc <decode+0x28c>
     864:	0f 1f 40 00          	nopl   0x0(%rax)
     868:	48 83 c4 58          	add    $0x58,%rsp
     86c:	b8 fd ff ff ff       	mov    $0xfffffffd,%eax
     871:	5b                   	pop    %rbx
     872:	5d                   	pop    %rbp
     873:	41 5c                	pop    %r12
     875:	41 5d                	pop    %r13
     877:	41 5e                	pop    %r14
     879:	41 5f                	pop    %r15
     87b:	c3                   	ret
     87c:	45 31 ff             	xor    %r15d,%r15d
     87f:	31 c0                	xor    %eax,%eax
     881:	bd 01 00 00 00       	mov    $0x1,%ebp
     886:	41 be 01 00 00 00    	mov    $0x1,%r14d
     88c:	48 c7 44 24 20 01 00 	movq   $0x1,0x20(%rsp)
     893:	00 00 
     895:	c7 44 24 44 00 00 00 	movl   $0x0,0x44(%rsp)
     89c:	00 
     89d:	e9 5a ff ff ff       	jmp    7fc <decode+0x2bc>
     8a2:	b8 ff ff ff ff       	mov    $0xffffffff,%eax
     8a7:	e9 85 fe ff ff       	jmp    731 <decode+0x1f1>
     8ac:	0f 1f 40 00          	nopl   0x0(%rax)

00000000000008b0 <ac_bits>:
     8b0:	f3 0f 1e fa          	endbr64
     8b4:	69 c7 8b 8a 01 00    	imul   $0x18a8b,%edi,%eax
     8ba:	05 3f 0d 03 00       	add    $0x30d3f,%eax
     8bf:	c1 e8 06             	shr    $0x6,%eax
     8c2:	48 69 c0 63 2d 3e 05 	imul   $0x53e2d63,%rax,%rax
     8c9:	48 c1 e8 26          	shr    $0x26,%rax
     8cd:	83 c0 03             	add    $0x3,%eax
     8d0:	c3                   	ret
     8d1:	66 66 2e 0f 1f 84 00 	data16 cs nopw 0x0(%rax,%rax,1)
     8d8:	00 00 00 00 
     8dc:	0f 1f 40 00          	nopl   0x0(%rax)

00000000000008e0 <ac_encode>:
     8e0:	f3 0f 1e fa          	endbr64
     8e4:	41 57                	push   %r15
     8e6:	41 56                	push   %r14
     8e8:	41 55                	push   %r13
     8ea:	41 54                	push   %r12
     8ec:	55                   	push   %rbp
     8ed:	53                   	push   %rbx
     8ee:	48 89 d3             	mov    %rdx,%rbx
     8f1:	8d 57 fc             	lea    -0x4(%rdi),%edx
     8f4:	48 83 ec 58          	sub    $0x58,%rsp
     8f8:	48 89 74 24 40       	mov    %rsi,0x40(%rsp)
     8fd:	83 fa 1a             	cmp    $0x1a,%edx
     900:	0f 87 fd 02 00 00    	ja     c03 <ac_encode+0x323>
     906:	48 85 f6             	test   %rsi,%rsi
     909:	0f 84 f4 02 00 00    	je     c03 <ac_encode+0x323>
     90f:	48 83 fe 01          	cmp    $0x1,%rsi
     913:	0f 84 df 02 00 00    	je     bf8 <ac_encode+0x318>
     919:	49 89 f6             	mov    %rsi,%r14
     91c:	31 c0                	xor    %eax,%eax
     91e:	66 90                	xchg   %ax,%ax
     920:	49 d1 ee             	shr    $1,%r14
     923:	83 c0 01             	add    $0x1,%eax
     926:	49 83 fe 01          	cmp    $0x1,%r14
     92a:	75 f4                	jne    920 <ac_encode+0x40>
     92c:	89 44 24 20          	mov    %eax,0x20(%rsp)
     930:	83 f8 14             	cmp    $0x14,%eax
     933:	0f 87 ca 02 00 00    	ja     c03 <ac_encode+0x323>
     939:	69 cf 8b 8a 01 00    	imul   $0x18a8b,%edi,%ecx
     93f:	4c 89 f2             	mov    %r14,%rdx
     942:	81 c1 3f 0d 03 00    	add    $0x30d3f,%ecx
     948:	c1 e9 06             	shr    $0x6,%ecx
     94b:	48 69 c9 63 2d 3e 05 	imul   $0x53e2d63,%rcx,%rcx
     952:	48 c1 e9 26          	shr    $0x26,%rcx
     956:	83 c1 03             	add    $0x3,%ecx
     959:	48 d3 e2             	shl    %cl,%rdx
     95c:	48 89 54 24 08       	mov    %rdx,0x8(%rsp)
     961:	85 c0                	test   %eax,%eax
     963:	0f 84 8f 02 00 00    	je     bf8 <ac_encode+0x318>
     969:	4c 89 f0             	mov    %r14,%rax
     96c:	89 f9                	mov    %edi,%ecx
     96e:	48 89 5c 24 48       	mov    %rbx,0x48(%rsp)
     973:	41 bf 01 00 00 00    	mov    $0x1,%r15d
     979:	48 c7 44 24 38 00 00 	movq   $0x0,0x38(%rsp)
     980:	00 00 
     982:	48 d3 e0             	shl    %cl,%rax
     985:	31 ed                	xor    %ebp,%ebp
     987:	48 89 44 24 30       	mov    %rax,0x30(%rsp)
     98c:	31 c0                	xor    %eax,%eax
     98e:	48 c7 44 24 28 00 00 	movq   $0x0,0x28(%rsp)
     995:	00 00 
     997:	48 89 c3             	mov    %rax,%rbx
     99a:	eb 3d                	jmp    9d9 <ac_encode+0xf9>
     99c:	0f 1f 40 00          	nopl   0x0(%rax)
     9a0:	4b 8d 04 36          	lea    (%r14,%r14,1),%rax
     9a4:	48 8d 4c 6d 00       	lea    0x0(%rbp,%rbp,2),%rcx
     9a9:	48 8d 14 28          	lea    (%rax,%rbp,1),%rdx
     9ad:	4c 01 f0             	add    %r14,%rax
     9b0:	4b 8d 34 7f          	lea    (%r15,%r15,2),%rsi
     9b4:	48 8d 2c 48          	lea    (%rax,%rcx,2),%rbp
     9b8:	48 8d 04 1b          	lea    (%rbx,%rbx,1),%rax
     9bc:	49 89 d6             	mov    %rdx,%r14
     9bf:	4a 8d 0c 38          	lea    (%rax,%r15,1),%rcx
     9c3:	48 01 d8             	add    %rbx,%rax
     9c6:	4c 8d 3c 70          	lea    (%rax,%rsi,2),%r15
     9ca:	8b 44 24 20          	mov    0x20(%rsp),%eax
     9ce:	48 89 cb             	mov    %rcx,%rbx
     9d1:	85 c0                	test   %eax,%eax
     9d3:	0f 84 e9 01 00 00    	je     bc2 <ac_encode+0x2e2>
     9d9:	83 6c 24 20 01       	subl   $0x1,0x20(%rsp)
     9de:	8b 44 24 20          	mov    0x20(%rsp),%eax
     9e2:	49 89 da             	mov    %rbx,%r10
     9e5:	45 31 db             	xor    %r11d,%r11d
     9e8:	48 8b 54 24 40       	mov    0x40(%rsp),%rdx
     9ed:	49 89 ec             	mov    %rbp,%r12
     9f0:	89 c1                	mov    %eax,%ecx
     9f2:	48 d3 ea             	shr    %cl,%rdx
     9f5:	48 89 d0             	mov    %rdx,%rax
     9f8:	83 e0 01             	and    $0x1,%eax
     9fb:	4d 01 f2             	add    %r14,%r10
     9fe:	49 83 d3 00          	adc    $0x0,%r11
     a02:	89 44 24 24          	mov    %eax,0x24(%rsp)
     a06:	45 31 ed             	xor    %r13d,%r13d
     a09:	4c 89 f0             	mov    %r14,%rax
     a0c:	4d 01 fc             	add    %r15,%r12
     a0f:	49 83 d5 00          	adc    $0x0,%r13
     a13:	49 f7 e7             	mul    %r15
     a16:	49 89 c0             	mov    %rax,%r8
     a19:	48 89 e8             	mov    %rbp,%rax
     a1c:	49 89 d1             	mov    %rdx,%r9
     a1f:	48 f7 e3             	mul    %rbx
     a22:	48 89 c6             	mov    %rax,%rsi
     a25:	4c 89 c0             	mov    %r8,%rax
     a28:	48 89 d7             	mov    %rdx,%rdi
     a2b:	4c 89 ca             	mov    %r9,%rdx
     a2e:	48 29 f0             	sub    %rsi,%rax
     a31:	48 8b 74 24 30       	mov    0x30(%rsp),%rsi
     a36:	48 19 fa             	sbb    %rdi,%rdx
     a39:	48 8b 7c 24 38       	mov    0x38(%rsp),%rdi
     a3e:	48 89 f1             	mov    %rsi,%rcx
     a41:	48 0f af ca          	imul   %rdx,%rcx
     a45:	48 89 fa             	mov    %rdi,%rdx
     a48:	48 0f af d0          	imul   %rax,%rdx
     a4c:	48 01 d1             	add    %rdx,%rcx
     a4f:	48 f7 e6             	mul    %rsi
     a52:	48 01 ca             	add    %rcx,%rdx
     a55:	48 89 c6             	mov    %rax,%rsi
     a58:	4c 89 d9             	mov    %r11,%rcx
     a5b:	48 89 d7             	mov    %rdx,%rdi
     a5e:	48 c1 e6 05          	shl    $0x5,%rsi
     a62:	48 0f a4 c7 05       	shld   $0x5,%rax,%rdi
     a67:	48 01 c6             	add    %rax,%rsi
     a6a:	4c 89 e0             	mov    %r12,%rax
     a6d:	48 11 d7             	adc    %rdx,%rdi
     a70:	4c 89 d2             	mov    %r10,%rdx
     a73:	4c 0f a4 d1 01       	shld   $0x1,%r10,%rcx
     a78:	48 0f a4 f7 04       	shld   $0x4,%rsi,%rdi
     a7d:	4c 01 d2             	add    %r10,%rdx
     a80:	48 c1 e6 04          	shl    $0x4,%rsi
     a84:	48 89 74 24 10       	mov    %rsi,0x10(%rsp)
     a89:	48 89 d6             	mov    %rdx,%rsi
     a8c:	4c 89 ea             	mov    %r13,%rdx
     a8f:	4c 01 d6             	add    %r10,%rsi
     a92:	48 89 7c 24 18       	mov    %rdi,0x18(%rsp)
     a97:	48 89 cf             	mov    %rcx,%rdi
     a9a:	4c 11 df             	adc    %r11,%rdi
     a9d:	48 c1 e0 02          	shl    $0x2,%rax
     aa1:	49 89 f0             	mov    %rsi,%r8
     aa4:	4c 89 e6             	mov    %r12,%rsi
     aa7:	4c 0f a4 e2 02       	shld   $0x2,%r12,%rdx
     aac:	4c 01 e0             	add    %r12,%rax
     aaf:	49 89 f9             	mov    %rdi,%r9
     ab2:	4c 89 ef             	mov    %r13,%rdi
     ab5:	4c 11 ea             	adc    %r13,%rdx
     ab8:	49 01 c0             	add    %rax,%r8
     abb:	4c 89 d0             	mov    %r10,%rax
     abe:	49 11 d1             	adc    %rdx,%r9
     ac1:	4c 89 da             	mov    %r11,%rdx
     ac4:	48 c1 e0 02          	shl    $0x2,%rax
     ac8:	4c 0f a4 d2 02       	shld   $0x2,%r10,%rdx
     acd:	4c 01 d0             	add    %r10,%rax
     ad0:	4c 89 c9             	mov    %r9,%rcx
     ad3:	4c 11 da             	adc    %r11,%rdx
     ad6:	4c 01 e6             	add    %r12,%rsi
     ad9:	4c 0f a4 e7 01       	shld   $0x1,%r12,%rdi
     ade:	4c 01 e6             	add    %r12,%rsi
     ae1:	4c 11 ef             	adc    %r13,%rdi
     ae4:	48 01 f0             	add    %rsi,%rax
     ae7:	48 11 fa             	adc    %rdi,%rdx
     aea:	48 0f af c8          	imul   %rax,%rcx
     aee:	49 0f af d0          	imul   %r8,%rdx
     af2:	48 01 d1             	add    %rdx,%rcx
     af5:	49 f7 e0             	mul    %r8
     af8:	48 01 ca             	add    %rcx,%rdx
     afb:	48 89 c6             	mov    %rax,%rsi
     afe:	48 89 d7             	mov    %rdx,%rdi
     b01:	48 01 c6             	add    %rax,%rsi
     b04:	48 0f a4 c7 01       	shld   $0x1,%rax,%rdi
     b09:	48 01 c6             	add    %rax,%rsi
     b0c:	48 11 d7             	adc    %rdx,%rdi
     b0f:	48 0f a4 f7 04       	shld   $0x4,%rsi,%rdi
     b14:	48 c1 e6 04          	shl    $0x4,%rsi
     b18:	48 01 c6             	add    %rax,%rsi
     b1b:	48 8b 44 24 10       	mov    0x10(%rsp),%rax
     b20:	48 11 d7             	adc    %rdx,%rdi
     b23:	48 8b 54 24 18       	mov    0x18(%rsp),%rdx
     b28:	48 39 c6             	cmp    %rax,%rsi
     b2b:	48 19 d7             	sbb    %rdx,%rdi
     b2e:	0f 83 ac 00 00 00    	jae    be0 <ac_encode+0x300>
     b34:	4c 8b 64 24 08       	mov    0x8(%rsp),%r12
     b39:	49 8d 34 2f          	lea    (%r15,%rbp,1),%rsi
     b3d:	4a 8d 3c 33          	lea    (%rbx,%r14,1),%rdi
     b41:	e8 00 00 00 00       	call   b46 <ac_encode+0x266>
			b42: R_X86_64_PLT32	ac_branch-0x4
     b46:	49 83 ec 01          	sub    $0x1,%r12
     b4a:	89 c0                	mov    %eax,%eax
     b4c:	49 0f af c4          	imul   %r12,%rax
     b50:	48 89 c1             	mov    %rax,%rcx
     b53:	48 c1 e9 10          	shr    $0x10,%rcx
     b57:	48 89 4c 24 08       	mov    %rcx,0x8(%rsp)
     b5c:	4c 39 e1             	cmp    %r12,%rcx
     b5f:	73 7f                	jae    be0 <ac_encode+0x300>
     b61:	48 3d ff ff 00 00    	cmp    $0xffff,%rax
     b67:	76 77                	jbe    be0 <ac_encode+0x300>
     b69:	8b 54 24 24          	mov    0x24(%rsp),%edx
     b6d:	48 83 44 24 28 01    	addq   $0x1,0x28(%rsp)
     b73:	85 d2                	test   %edx,%edx
     b75:	0f 84 25 fe ff ff    	je     9a0 <ac_encode+0xc0>
     b7b:	48 8b 54 24 08       	mov    0x8(%rsp),%rdx
     b80:	48 8d 44 2d 00       	lea    0x0(%rbp,%rbp,1),%rax
     b85:	48 01 54 24 28       	add    %rdx,0x28(%rsp)
     b8a:	4b 8d 0c 76          	lea    (%r14,%r14,2),%rcx
     b8e:	48 8d 3c 5b          	lea    (%rbx,%rbx,2),%rdi
     b92:	49 29 d4             	sub    %rdx,%r12
     b95:	48 8d 14 28          	lea    (%rax,%rbp,1),%rdx
     b99:	4a 8d 2c 30          	lea    (%rax,%r14,1),%rbp
     b9d:	4b 8d 04 3f          	lea    (%r15,%r15,1),%rax
     ba1:	4c 89 64 24 08       	mov    %r12,0x8(%rsp)
     ba6:	4c 8d 34 4a          	lea    (%rdx,%rcx,2),%r14
     baa:	4a 8d 34 38          	lea    (%rax,%r15,1),%rsi
     bae:	4c 8d 3c 18          	lea    (%rax,%rbx,1),%r15
     bb2:	8b 44 24 20          	mov    0x20(%rsp),%eax
     bb6:	48 8d 1c 7e          	lea    (%rsi,%rdi,2),%rbx
     bba:	85 c0                	test   %eax,%eax
     bbc:	0f 85 17 fe ff ff    	jne    9d9 <ac_encode+0xf9>
     bc2:	48 8b 5c 24 48       	mov    0x48(%rsp),%rbx
     bc7:	48 8b 44 24 28       	mov    0x28(%rsp),%rax
     bcc:	48 89 03             	mov    %rax,(%rbx)
     bcf:	48 83 c4 58          	add    $0x58,%rsp
     bd3:	31 c0                	xor    %eax,%eax
     bd5:	5b                   	pop    %rbx
     bd6:	5d                   	pop    %rbp
     bd7:	41 5c                	pop    %r12
     bd9:	41 5d                	pop    %r13
     bdb:	41 5e                	pop    %r14
     bdd:	41 5f                	pop    %r15
     bdf:	c3                   	ret
     be0:	b8 fe ff ff ff       	mov    $0xfffffffe,%eax
     be5:	48 83 c4 58          	add    $0x58,%rsp
     be9:	5b                   	pop    %rbx
     bea:	5d                   	pop    %rbp
     beb:	41 5c                	pop    %r12
     bed:	41 5d                	pop    %r13
     bef:	41 5e                	pop    %r14
     bf1:	41 5f                	pop    %r15
     bf3:	c3                   	ret
     bf4:	0f 1f 40 00          	nopl   0x0(%rax)
     bf8:	48 c7 44 24 28 00 00 	movq   $0x0,0x28(%rsp)
     bff:	00 00 
     c01:	eb c4                	jmp    bc7 <ac_encode+0x2e7>
     c03:	b8 ff ff ff ff       	mov    $0xffffffff,%eax
     c08:	eb db                	jmp    be5 <ac_encode+0x305>
     c0a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)

0000000000000c10 <ac_decode>:
     c10:	f3 0f 1e fa          	endbr64
     c14:	53                   	push   %rbx
     c15:	48 89 d3             	mov    %rdx,%rbx
     c18:	48 83 ec 40          	sub    $0x40,%rsp
     c1c:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
     c23:	00 00 
     c25:	48 89 44 24 38       	mov    %rax,0x38(%rsp)
     c2a:	31 c0                	xor    %eax,%eax
     c2c:	48 89 e2             	mov    %rsp,%rdx
     c2f:	e8 0c f9 ff ff       	call   540 <decode>
     c34:	85 c0                	test   %eax,%eax
     c36:	75 08                	jne    c40 <ac_decode+0x30>
     c38:	48 8b 54 24 28       	mov    0x28(%rsp),%rdx
     c3d:	48 89 13             	mov    %rdx,(%rbx)
     c40:	48 8b 54 24 38       	mov    0x38(%rsp),%rdx
     c45:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
     c4c:	00 00 
     c4e:	75 06                	jne    c56 <ac_decode+0x46>
     c50:	48 83 c4 40          	add    $0x40,%rsp
     c54:	5b                   	pop    %rbx
     c55:	c3                   	ret
     c56:	e8 00 00 00 00       	call   c5b <ac_decode+0x4b>
			c57: R_X86_64_PLT32	__stack_chk_fail-0x4
     c5b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000000c60 <ac_output>:
     c60:	f3 0f 1e fa          	endbr64
     c64:	41 54                	push   %r12
     c66:	55                   	push   %rbp
     c67:	89 fd                	mov    %edi,%ebp
     c69:	53                   	push   %rbx
     c6a:	48 83 ec 40          	sub    $0x40,%rsp
     c6e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
     c75:	00 00 
     c77:	48 89 44 24 38       	mov    %rax,0x38(%rsp)
     c7c:	31 c0                	xor    %eax,%eax
     c7e:	8d 47 fc             	lea    -0x4(%rdi),%eax
     c81:	83 f8 1a             	cmp    $0x1a,%eax
     c84:	0f 87 ab 00 00 00    	ja     d35 <ac_output+0xd5>
     c8a:	41 89 d0             	mov    %edx,%r8d
     c8d:	41 8d 40 ff          	lea    -0x1(%r8),%eax
     c91:	83 f8 03             	cmp    $0x3,%eax
     c94:	0f 87 9b 00 00 00    	ja     d35 <ac_output+0xd5>
     c9a:	69 c7 8b 8a 01 00    	imul   $0x18a8b,%edi,%eax
     ca0:	42 8d 14 c5 00 00 00 	lea    0x0(,%r8,8),%edx
     ca7:	00 
     ca8:	05 3f 0d 03 00       	add    $0x30d3f,%eax
     cad:	c1 e8 06             	shr    $0x6,%eax
     cb0:	48 69 c0 63 2d 3e 05 	imul   $0x53e2d63,%rax,%rax
     cb7:	48 c1 e8 26          	shr    $0x26,%rax
     cbb:	83 c0 03             	add    $0x3,%eax
     cbe:	39 c2                	cmp    %eax,%edx
     cc0:	72 73                	jb     d35 <ac_output+0xd5>
     cc2:	48 89 f7             	mov    %rsi,%rdi
     cc5:	48 89 cb             	mov    %rcx,%rbx
     cc8:	31 c0                	xor    %eax,%eax
     cca:	31 f6                	xor    %esi,%esi
     ccc:	0f 1f 40 00          	nopl   0x0(%rax)
     cd0:	0f b6 14 07          	movzbl (%rdi,%rax,1),%edx
     cd4:	8d 0c c5 00 00 00 00 	lea    0x0(,%rax,8),%ecx
     cdb:	48 83 c0 01          	add    $0x1,%rax
     cdf:	48 d3 e2             	shl    %cl,%rdx
     ce2:	48 09 d6             	or     %rdx,%rsi
     ce5:	49 39 c0             	cmp    %rax,%r8
     ce8:	75 e6                	jne    cd0 <ac_output+0x70>
     cea:	48 89 e2             	mov    %rsp,%rdx
     ced:	89 ef                	mov    %ebp,%edi
     cef:	e8 4c f8 ff ff       	call   540 <decode>
     cf4:	41 89 c4             	mov    %eax,%r12d
     cf7:	85 c0                	test   %eax,%eax
     cf9:	75 42                	jne    d3d <ac_output+0xdd>
     cfb:	48 8b 4c 24 18       	mov    0x18(%rsp),%rcx
     d00:	48 8b 54 24 10       	mov    0x10(%rsp),%rdx
     d05:	49 89 d9             	mov    %rbx,%r9
     d08:	41 89 e8             	mov    %ebp,%r8d
     d0b:	48 8b 74 24 08       	mov    0x8(%rsp),%rsi
     d10:	48 8b 3c 24          	mov    (%rsp),%rdi
     d14:	e8 b7 f4 ff ff       	call   1d0 <output.isra.0>
     d19:	48 8b 44 24 38       	mov    0x38(%rsp),%rax
     d1e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
     d25:	00 00 
     d27:	75 1c                	jne    d45 <ac_output+0xe5>
     d29:	48 83 c4 40          	add    $0x40,%rsp
     d2d:	44 89 e0             	mov    %r12d,%eax
     d30:	5b                   	pop    %rbx
     d31:	5d                   	pop    %rbp
     d32:	41 5c                	pop    %r12
     d34:	c3                   	ret
     d35:	41 bc ff ff ff ff    	mov    $0xffffffff,%r12d
     d3b:	eb dc                	jmp    d19 <ac_output+0xb9>
     d3d:	41 bc fe ff ff ff    	mov    $0xfffffffe,%r12d
     d43:	eb d4                	jmp    d19 <ac_output+0xb9>
     d45:	e8 00 00 00 00       	call   d4a <ac_output+0xea>
			d46: R_X86_64_PLT32	__stack_chk_fail-0x4
     d4a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)

0000000000000d50 <ac_advance>:
     d50:	f3 0f 1e fa          	endbr64
     d54:	41 57                	push   %r15
     d56:	89 d0                	mov    %edx,%eax
     d58:	41 56                	push   %r14
     d5a:	41 55                	push   %r13
     d5c:	41 54                	push   %r12
     d5e:	55                   	push   %rbp
     d5f:	53                   	push   %rbx
     d60:	48 81 ec a8 00 00 00 	sub    $0xa8,%rsp
     d67:	4c 89 44 24 40       	mov    %r8,0x40(%rsp)
     d6c:	64 48 8b 14 25 28 00 	mov    %fs:0x28,%rdx
     d73:	00 00 
     d75:	48 89 94 24 98 00 00 	mov    %rdx,0x98(%rsp)
     d7c:	00 
     d7d:	31 d2                	xor    %edx,%edx
     d7f:	8d 57 fc             	lea    -0x4(%rdi),%edx
     d82:	83 fa 1a             	cmp    $0x1a,%edx
     d85:	0f 87 1f 04 00 00    	ja     11aa <ac_advance+0x45a>
     d8b:	8d 50 ff             	lea    -0x1(%rax),%edx
     d8e:	83 fa 03             	cmp    $0x3,%edx
     d91:	0f 87 13 04 00 00    	ja     11aa <ac_advance+0x45a>
     d97:	44 69 f7 8b 8a 01 00 	imul   $0x18a8b,%edi,%r14d
     d9e:	8d 14 c5 00 00 00 00 	lea    0x0(,%rax,8),%edx
     da5:	89 fb                	mov    %edi,%ebx
     da7:	41 81 c6 3f 0d 03 00 	add    $0x30d3f,%r14d
     dae:	41 c1 ee 06          	shr    $0x6,%r14d
     db2:	4d 69 f6 63 2d 3e 05 	imul   $0x53e2d63,%r14,%r14
     db9:	49 c1 ee 26          	shr    $0x26,%r14
     dbd:	41 83 c6 03          	add    $0x3,%r14d
     dc1:	44 39 f2             	cmp    %r14d,%edx
     dc4:	0f 82 e0 03 00 00    	jb     11aa <ac_advance+0x45a>
     dca:	89 cd                	mov    %ecx,%ebp
     dcc:	83 f9 01             	cmp    $0x1,%ecx
     dcf:	0f 87 d5 03 00 00    	ja     11aa <ac_advance+0x45a>
     dd5:	49 89 f5             	mov    %rsi,%r13
     dd8:	41 89 c4             	mov    %eax,%r12d
     ddb:	31 f6                	xor    %esi,%esi
     ddd:	31 c0                	xor    %eax,%eax
     ddf:	90                   	nop
     de0:	41 0f b6 54 05 00    	movzbl 0x0(%r13,%rax,1),%edx
     de6:	8d 0c c5 00 00 00 00 	lea    0x0(,%rax,8),%ecx
     ded:	48 83 c0 01          	add    $0x1,%rax
     df1:	48 d3 e2             	shl    %cl,%rdx
     df4:	48 09 d6             	or     %rdx,%rsi
     df7:	4c 39 e0             	cmp    %r12,%rax
     dfa:	75 e4                	jne    de0 <ac_advance+0x90>
     dfc:	48 8d 54 24 60       	lea    0x60(%rsp),%rdx
     e01:	89 df                	mov    %ebx,%edi
     e03:	e8 38 f7 ff ff       	call   540 <decode>
     e08:	89 44 24 50          	mov    %eax,0x50(%rsp)
     e0c:	85 c0                	test   %eax,%eax
     e0e:	0f 85 a0 03 00 00    	jne    11b4 <ac_advance+0x464>
     e14:	44 89 f1             	mov    %r14d,%ecx
     e17:	b8 01 00 00 00       	mov    $0x1,%eax
     e1c:	be 01 00 00 00       	mov    $0x1,%esi
     e21:	bf 01 00 00 00       	mov    $0x1,%edi
     e26:	48 d3 e0             	shl    %cl,%rax
     e29:	4c 8d 70 ff          	lea    -0x1(%rax),%r14
     e2d:	e8 00 00 00 00       	call   e32 <ac_advance+0xe2>
			e2e: R_X86_64_PLT32	ac_branch-0x4
     e32:	89 c0                	mov    %eax,%eax
     e34:	49 0f af c6          	imul   %r14,%rax
     e38:	48 89 c7             	mov    %rax,%rdi
     e3b:	48 c1 ef 10          	shr    $0x10,%rdi
     e3f:	48 89 fe             	mov    %rdi,%rsi
     e42:	48 3d ff ff 00 00    	cmp    $0xffff,%rax
     e48:	0f 86 52 03 00 00    	jbe    11a0 <ac_advance+0x450>
     e4e:	4c 39 f7             	cmp    %r14,%rdi
     e51:	0f 83 49 03 00 00    	jae    11a0 <ac_advance+0x450>
     e57:	85 ed                	test   %ebp,%ebp
     e59:	0f 85 11 03 00 00    	jne    1170 <ac_advance+0x420>
     e5f:	41 bf 06 00 00 00    	mov    $0x6,%r15d
     e65:	bf 01 00 00 00       	mov    $0x1,%edi
     e6a:	bd 02 00 00 00       	mov    $0x2,%ebp
     e6f:	48 c7 44 24 30 01 00 	movq   $0x1,0x30(%rsp)
     e76:	00 00 
     e78:	41 be 03 00 00 00    	mov    $0x3,%r14d
     e7e:	8b 84 24 80 00 00 00 	mov    0x80(%rsp),%eax
     e85:	89 44 24 0c          	mov    %eax,0xc(%rsp)
     e89:	85 c0                	test   %eax,%eax
     e8b:	0f 84 72 02 00 00    	je     1103 <ac_advance+0x3b3>
     e91:	b8 01 00 00 00       	mov    $0x1,%eax
     e96:	89 d9                	mov    %ebx,%ecx
     e98:	4c 89 64 24 48       	mov    %r12,0x48(%rsp)
     e9d:	49 89 fc             	mov    %rdi,%r12
     ea0:	48 d3 e0             	shl    %cl,%rax
     ea3:	4c 89 6c 24 58       	mov    %r13,0x58(%rsp)
     ea8:	49 89 f5             	mov    %rsi,%r13
     eab:	48 89 44 24 20       	mov    %rax,0x20(%rsp)
     eb0:	48 8b 84 24 88 00 00 	mov    0x88(%rsp),%rax
     eb7:	00 
     eb8:	48 c7 44 24 28 00 00 	movq   $0x0,0x28(%rsp)
     ebf:	00 00 
     ec1:	48 89 44 24 38       	mov    %rax,0x38(%rsp)
     ec6:	89 5c 24 54          	mov    %ebx,0x54(%rsp)
     eca:	eb 49                	jmp    f15 <ac_advance+0x1c5>
     ecc:	0f 1f 40 00          	nopl   0x0(%rax)
     ed0:	4c 01 e8             	add    %r13,%rax
     ed3:	4c 29 ea             	sub    %r13,%rdx
     ed6:	44 8b 44 24 0c       	mov    0xc(%rsp),%r8d
     edb:	48 8d 4c 6d 00       	lea    0x0(%rbp,%rbp,2),%rcx
     ee0:	48 89 44 24 30       	mov    %rax,0x30(%rsp)
     ee5:	4b 8d 04 36          	lea    (%r14,%r14,1),%rax
     ee9:	49 89 d5             	mov    %rdx,%r13
     eec:	4b 8d 3c 64          	lea    (%r12,%r12,2),%rdi
     ef0:	4a 8d 14 30          	lea    (%rax,%r14,1),%rdx
     ef4:	4c 8d 34 28          	lea    (%rax,%rbp,1),%r14
     ef8:	4b 8d 04 3f          	lea    (%r15,%r15,1),%rax
     efc:	48 8d 2c 4a          	lea    (%rdx,%rcx,2),%rbp
     f00:	4a 8d 34 38          	lea    (%rax,%r15,1),%rsi
     f04:	4e 8d 3c 20          	lea    (%rax,%r12,1),%r15
     f08:	4c 8d 24 7e          	lea    (%rsi,%rdi,2),%r12
     f0c:	45 85 c0             	test   %r8d,%r8d
     f0f:	0f 84 dd 01 00 00    	je     10f2 <ac_advance+0x3a2>
     f15:	4c 89 e6             	mov    %r12,%rsi
     f18:	31 ff                	xor    %edi,%edi
     f1a:	4d 89 f2             	mov    %r14,%r10
     f1d:	48 89 e8             	mov    %rbp,%rax
     f20:	48 01 ee             	add    %rbp,%rsi
     f23:	48 83 d7 00          	adc    $0x0,%rdi
     f27:	45 31 db             	xor    %r11d,%r11d
     f2a:	4d 01 fa             	add    %r15,%r10
     f2d:	49 83 d3 00          	adc    $0x0,%r11
     f31:	49 f7 e7             	mul    %r15
     f34:	48 89 c1             	mov    %rax,%rcx
     f37:	4c 89 f0             	mov    %r14,%rax
     f3a:	48 89 d3             	mov    %rdx,%rbx
     f3d:	49 f7 e4             	mul    %r12
     f40:	49 89 c0             	mov    %rax,%r8
     f43:	48 89 c8             	mov    %rcx,%rax
     f46:	49 89 d1             	mov    %rdx,%r9
     f49:	48 89 da             	mov    %rbx,%rdx
     f4c:	4c 29 c0             	sub    %r8,%rax
     f4f:	4c 8b 44 24 20       	mov    0x20(%rsp),%r8
     f54:	4c 19 ca             	sbb    %r9,%rdx
     f57:	4c 8b 4c 24 28       	mov    0x28(%rsp),%r9
     f5c:	48 89 44 24 10       	mov    %rax,0x10(%rsp)
     f61:	4c 89 c1             	mov    %r8,%rcx
     f64:	48 8b 5c 24 10       	mov    0x10(%rsp),%rbx
     f69:	48 89 54 24 18       	mov    %rdx,0x18(%rsp)
     f6e:	48 0f af ca          	imul   %rdx,%rcx
     f72:	4c 89 ca             	mov    %r9,%rdx
     f75:	48 0f af d3          	imul   %rbx,%rdx
     f79:	48 89 d8             	mov    %rbx,%rax
     f7c:	48 01 d1             	add    %rdx,%rcx
     f7f:	49 f7 e0             	mul    %r8
     f82:	48 01 ca             	add    %rcx,%rdx
     f85:	48 89 c1             	mov    %rax,%rcx
     f88:	48 89 d3             	mov    %rdx,%rbx
     f8b:	48 c1 e1 05          	shl    $0x5,%rcx
     f8f:	48 0f a4 c3 05       	shld   $0x5,%rax,%rbx
     f94:	48 01 c1             	add    %rax,%rcx
     f97:	48 89 f0             	mov    %rsi,%rax
     f9a:	48 11 d3             	adc    %rdx,%rbx
     f9d:	48 89 fa             	mov    %rdi,%rdx
     fa0:	48 01 f0             	add    %rsi,%rax
     fa3:	48 0f a4 cb 04       	shld   $0x4,%rcx,%rbx
     fa8:	48 c1 e1 04          	shl    $0x4,%rcx
     fac:	48 0f a4 f2 01       	shld   $0x1,%rsi,%rdx
     fb1:	49 89 c8             	mov    %rcx,%r8
     fb4:	48 01 f0             	add    %rsi,%rax
     fb7:	49 89 d9             	mov    %rbx,%r9
     fba:	4c 89 d1             	mov    %r10,%rcx
     fbd:	4c 89 db             	mov    %r11,%rbx
     fc0:	48 11 fa             	adc    %rdi,%rdx
     fc3:	48 0f a4 cb 02       	shld   $0x2,%rcx,%rbx
     fc8:	48 c1 e1 02          	shl    $0x2,%rcx
     fcc:	4c 01 d1             	add    %r10,%rcx
     fcf:	4c 11 db             	adc    %r11,%rbx
     fd2:	48 01 c8             	add    %rcx,%rax
     fd5:	48 89 f1             	mov    %rsi,%rcx
     fd8:	48 11 da             	adc    %rbx,%rdx
     fdb:	48 89 fb             	mov    %rdi,%rbx
     fde:	48 0f a4 cb 02       	shld   $0x2,%rcx,%rbx
     fe3:	48 c1 e1 02          	shl    $0x2,%rcx
     fe7:	48 01 f1             	add    %rsi,%rcx
     fea:	4c 89 d6             	mov    %r10,%rsi
     fed:	48 11 fb             	adc    %rdi,%rbx
     ff0:	4c 89 df             	mov    %r11,%rdi
     ff3:	48 0f a4 f7 01       	shld   $0x1,%rsi,%rdi
     ff8:	48 01 f6             	add    %rsi,%rsi
     ffb:	4c 01 d6             	add    %r10,%rsi
     ffe:	4c 11 df             	adc    %r11,%rdi
    1001:	48 01 f1             	add    %rsi,%rcx
    1004:	48 89 d6             	mov    %rdx,%rsi
    1007:	48 11 fb             	adc    %rdi,%rbx
    100a:	48 0f af f1          	imul   %rcx,%rsi
    100e:	48 89 da             	mov    %rbx,%rdx
    1011:	48 0f af d0          	imul   %rax,%rdx
    1015:	48 01 d6             	add    %rdx,%rsi
    1018:	48 f7 e1             	mul    %rcx
    101b:	48 01 f2             	add    %rsi,%rdx
    101e:	48 89 c1             	mov    %rax,%rcx
    1021:	48 01 c1             	add    %rax,%rcx
    1024:	48 89 d3             	mov    %rdx,%rbx
    1027:	48 0f a4 c3 01       	shld   $0x1,%rax,%rbx
    102c:	48 01 c1             	add    %rax,%rcx
    102f:	48 11 d3             	adc    %rdx,%rbx
    1032:	48 0f a4 cb 04       	shld   $0x4,%rcx,%rbx
    1037:	48 c1 e1 04          	shl    $0x4,%rcx
    103b:	48 01 c1             	add    %rax,%rcx
    103e:	48 11 d3             	adc    %rdx,%rbx
    1041:	4c 39 c1             	cmp    %r8,%rcx
    1044:	48 89 d8             	mov    %rbx,%rax
    1047:	4c 19 c8             	sbb    %r9,%rax
    104a:	0f 83 a2 00 00 00    	jae    10f2 <ac_advance+0x3a2>
    1050:	83 6c 24 0c 01       	subl   $0x1,0xc(%rsp)
    1055:	8b 44 24 0c          	mov    0xc(%rsp),%eax
    1059:	49 8d 55 ff          	lea    -0x1(%r13),%rdx
    105d:	4b 8d 34 37          	lea    (%r15,%r14,1),%rsi
    1061:	48 8b 5c 24 38       	mov    0x38(%rsp),%rbx
    1066:	4a 8d 7c 25 00       	lea    0x0(%rbp,%r12,1),%rdi
    106b:	48 89 54 24 10       	mov    %rdx,0x10(%rsp)
    1070:	89 c1                	mov    %eax,%ecx
    1072:	48 d3 eb             	shr    %cl,%rbx
    1075:	e8 00 00 00 00       	call   107a <ac_advance+0x32a>
			1076: R_X86_64_PLT32	ac_branch-0x4
    107a:	48 8b 54 24 10       	mov    0x10(%rsp),%rdx
    107f:	83 e3 01             	and    $0x1,%ebx
    1082:	89 c0                	mov    %eax,%eax
    1084:	48 0f af c2          	imul   %rdx,%rax
    1088:	49 89 c5             	mov    %rax,%r13
    108b:	49 c1 ed 10          	shr    $0x10,%r13
    108f:	48 3d ff ff 00 00    	cmp    $0xffff,%rax
    1095:	0f 86 05 01 00 00    	jbe    11a0 <ac_advance+0x450>
    109b:	49 39 d5             	cmp    %rdx,%r13
    109e:	0f 83 fc 00 00 00    	jae    11a0 <ac_advance+0x450>
    10a4:	48 83 44 24 30 01    	addq   $0x1,0x30(%rsp)
    10aa:	48 8b 44 24 30       	mov    0x30(%rsp),%rax
    10af:	85 db                	test   %ebx,%ebx
    10b1:	0f 85 19 fe ff ff    	jne    ed0 <ac_advance+0x180>
    10b7:	48 8d 44 2d 00       	lea    0x0(%rbp,%rbp,1),%rax
    10bc:	4b 8d 0c 76          	lea    (%r14,%r14,2),%rcx
    10c0:	44 8b 44 24 0c       	mov    0xc(%rsp),%r8d
    10c5:	4a 8d 14 30          	lea    (%rax,%r14,1),%rdx
    10c9:	48 01 e8             	add    %rbp,%rax
    10cc:	4b 8d 34 7f          	lea    (%r15,%r15,2),%rsi
    10d0:	4c 8d 34 48          	lea    (%rax,%rcx,2),%r14
    10d4:	4b 8d 04 24          	lea    (%r12,%r12,1),%rax
    10d8:	48 89 d5             	mov    %rdx,%rbp
    10db:	4a 8d 0c 38          	lea    (%rax,%r15,1),%rcx
    10df:	4c 01 e0             	add    %r12,%rax
    10e2:	4c 8d 3c 70          	lea    (%rax,%rsi,2),%r15
    10e6:	49 89 cc             	mov    %rcx,%r12
    10e9:	45 85 c0             	test   %r8d,%r8d
    10ec:	0f 85 23 fe ff ff    	jne    f15 <ac_advance+0x1c5>
    10f2:	4c 89 e7             	mov    %r12,%rdi
    10f5:	8b 5c 24 54          	mov    0x54(%rsp),%ebx
    10f9:	4c 8b 64 24 48       	mov    0x48(%rsp),%r12
    10fe:	4c 8b 6c 24 58       	mov    0x58(%rsp),%r13
    1103:	48 8b 74 24 30       	mov    0x30(%rsp),%rsi
    1108:	31 c0                	xor    %eax,%eax
    110a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
    1110:	8d 0c c5 00 00 00 00 	lea    0x0(,%rax,8),%ecx
    1117:	48 89 f2             	mov    %rsi,%rdx
    111a:	48 d3 ea             	shr    %cl,%rdx
    111d:	41 88 54 05 00       	mov    %dl,0x0(%r13,%rax,1)
    1122:	48 83 c0 01          	add    $0x1,%rax
    1126:	4c 39 e0             	cmp    %r12,%rax
    1129:	75 e5                	jne    1110 <ac_advance+0x3c0>
    112b:	4c 8b 4c 24 40       	mov    0x40(%rsp),%r9
    1130:	48 89 fa             	mov    %rdi,%rdx
    1133:	41 89 d8             	mov    %ebx,%r8d
    1136:	4c 89 f9             	mov    %r15,%rcx
    1139:	4c 89 f6             	mov    %r14,%rsi
    113c:	48 89 ef             	mov    %rbp,%rdi
    113f:	e8 8c f0 ff ff       	call   1d0 <output.isra.0>
    1144:	48 8b 84 24 98 00 00 	mov    0x98(%rsp),%rax
    114b:	00 
    114c:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    1153:	00 00 
    1155:	75 67                	jne    11be <ac_advance+0x46e>
    1157:	8b 44 24 50          	mov    0x50(%rsp),%eax
    115b:	48 81 c4 a8 00 00 00 	add    $0xa8,%rsp
    1162:	5b                   	pop    %rbx
    1163:	5d                   	pop    %rbp
    1164:	41 5c                	pop    %r12
    1166:	41 5d                	pop    %r13
    1168:	41 5e                	pop    %r14
    116a:	41 5f                	pop    %r15
    116c:	c3                   	ret
    116d:	0f 1f 00             	nopl   (%rax)
    1170:	48 8d 7f 01          	lea    0x1(%rdi),%rdi
    1174:	49 29 f6             	sub    %rsi,%r14
    1177:	41 bf 02 00 00 00    	mov    $0x2,%r15d
    117d:	bd 06 00 00 00       	mov    $0x6,%ebp
    1182:	48 89 7c 24 30       	mov    %rdi,0x30(%rsp)
    1187:	4c 89 f6             	mov    %r14,%rsi
    118a:	bf 03 00 00 00       	mov    $0x3,%edi
    118f:	41 be 01 00 00 00    	mov    $0x1,%r14d
    1195:	e9 e4 fc ff ff       	jmp    e7e <ac_advance+0x12e>
    119a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
    11a0:	c7 44 24 50 fd ff ff 	movl   $0xfffffffd,0x50(%rsp)
    11a7:	ff 
    11a8:	eb 9a                	jmp    1144 <ac_advance+0x3f4>
    11aa:	c7 44 24 50 ff ff ff 	movl   $0xffffffff,0x50(%rsp)
    11b1:	ff 
    11b2:	eb 90                	jmp    1144 <ac_advance+0x3f4>
    11b4:	c7 44 24 50 fe ff ff 	movl   $0xfffffffe,0x50(%rsp)
    11bb:	ff 
    11bc:	eb 86                	jmp    1144 <ac_advance+0x3f4>
    11be:	e8 00 00 00 00       	call   11c3 <ac_advance+0x473>
			11bf: R_X86_64_PLT32	__stack_chk_fail-0x4
