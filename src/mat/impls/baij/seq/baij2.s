	.file	"baij2.c"
	.section	.debug_abbrev,"",@progbits
.Ldebug_abbrev0:
	.section	.debug_info,"",@progbits
.Ldebug_info0:
	.section	.debug_line,"",@progbits
.Ldebug_line0:
	.text
.Ltext0:
	.type	PetscAbsScalar, @function
PetscAbsScalar:
.LFB0:
	.file 1 "/home/dpnkarthik/petsc-rnet/include/petscmath.h"
	.loc 1 134 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	vmovsd	%xmm0, -8(%rbp)
	.loc 1 134 0
	vxorpd	%xmm0, %xmm0, %xmm0
	vucomisd	-8(%rbp), %xmm0
	seta	%al
	testb	%al, %al
	je	.L2
	vmovsd	-8(%rbp), %xmm1
	vmovsd	.LC1(%rip), %xmm0
	vxorpd	%xmm1, %xmm0, %xmm0
	jmp	.L3
.L2:
	vmovsd	-8(%rbp), %xmm0
.L3:
	vmovsd	%xmm0, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, -16(%rbp)
	vmovsd	-16(%rbp), %xmm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE0:
	.size	PetscAbsScalar, .-PetscAbsScalar
	.section	.rodata
	.align 8
.LC2:
	.string	"/home/dpnkarthik/petsc-rnet/include/petscsys.h"
.LC3:
	.string	"src/mat/impls/baij/seq/"
	.text
	.type	PetscMemcpy, @function
PetscMemcpy:
.LFB6:
	.file 2 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
	.loc 2 1763 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	.loc 2 1771 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L11
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L12
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.12688, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC2, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1771, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L10
.L11:
	nop
	jmp	.L10
.L12:
	nop
.L10:
	.loc 2 1773 0
	movq	-8(%rbp), %rax
	cmpq	-16(%rbp), %rax
	je	.L7
	.loc 2 1798 0
	movq	-24(%rbp), %rdx
	movq	-16(%rbp), %rcx
	movq	-8(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	memcpy
.L7:
	.loc 2 1801 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L8
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L8
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L8:
	movl	$0, %eax
	.loc 2 1802 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE6:
	.size	PetscMemcpy, .-PetscMemcpy
	.type	PetscMemzero, @function
PetscMemzero:
.LFB7:
	.loc 2 1827 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	.loc 2 1828 0
	cmpq	$0, -16(%rbp)
	je	.L14
	.loc 2 1847 0
	movq	-16(%rbp), %rdx
	movq	-8(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	memset
.L14:
	.loc 2 1853 0
	movl	$0, %eax
	.loc 2 1854 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE7:
	.size	PetscMemzero, .-PetscMemzero
	.section	.rodata
	.align 8
.LC4:
	.string	"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
.LC5:
	.string	" "
	.text
	.type	VecGetArrayRead, @function
VecGetArrayRead:
.LFB15:
	.file 3 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
	.loc 3 278 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	.loc 3 281 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L26
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L27
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.15175, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC4, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$281, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L25
.L26:
	nop
	jmp	.L25
.L27:
	nop
.L25:
	.loc 3 282 0
	movq	-24(%rbp), %rax
	movl	816(%rax), %eax
	testl	%eax, %eax
	je	.L18
	.loc 3 284 0
	movq	-24(%rbp), %rax
	movl	820(%rax), %eax
	cmpl	$1, %eax
	je	.L19
	movq	-24(%rbp), %rax
	movq	464(%rax), %rax
	movq	(%rax), %rax
	testq	%rax, %rax
	jne	.L20
.L19:
	.loc 3 285 0
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	VecCUSPCopyFromGPU
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L20
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC4, %ecx
	movl	$__func__.15175, %edx
	movl	$285, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L21
.L20:
	.loc 3 288 0
	movq	-24(%rbp), %rax
	movq	464(%rax), %rax
	movq	(%rax), %rdx
	movq	-32(%rbp), %rax
	movq	%rdx, (%rax)
	jmp	.L22
.L18:
	.loc 3 290 0
	movq	-24(%rbp), %rax
	movq	448(%rax), %rax
	movq	184(%rax), %rcx
	movq	-32(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	*%rcx
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L22
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC4, %ecx
	movl	$__func__.15175, %edx
	movl	$290, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L21
.L22:
	.loc 3 292 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L23
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L23
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L23:
	movl	$0, %eax
.L21:
	.loc 3 293 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE15:
	.size	VecGetArrayRead, .-VecGetArrayRead
	.type	VecRestoreArrayRead, @function
VecRestoreArrayRead:
.LFB16:
	.loc 3 298 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	.loc 3 301 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L36
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L37
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.15257, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC4, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$301, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L35
.L36:
	nop
	jmp	.L35
.L37:
	nop
.L35:
	.loc 3 302 0
	movq	-24(%rbp), %rax
	movl	816(%rax), %eax
	testl	%eax, %eax
	je	.L30
	.loc 3 304 0
	movq	-24(%rbp), %rax
	movl	820(%rax), %eax
	testl	%eax, %eax
	je	.L31
	.loc 3 305 0
	movq	-24(%rbp), %rax
	movl	$3, 820(%rax)
	jmp	.L31
.L30:
	.loc 3 309 0
	movq	-24(%rbp), %rax
	movq	448(%rax), %rax
	movq	208(%rax), %rcx
	movq	-32(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	*%rcx
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L31
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC4, %ecx
	movl	$__func__.15257, %edx
	movl	$309, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L32
.L31:
	.loc 3 311 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L33
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L33
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L33:
	movl	$0, %eax
.L32:
	.loc 3 312 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE16:
	.size	VecRestoreArrayRead, .-VecRestoreArrayRead
	.type	VecGetArray, @function
VecGetArray:
.LFB17:
	.loc 3 317 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	.loc 3 320 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L48
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L49
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.15326, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC4, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$320, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L47
.L48:
	nop
	jmp	.L47
.L49:
	nop
.L47:
	.loc 3 321 0
	movq	-24(%rbp), %rax
	movl	816(%rax), %eax
	testl	%eax, %eax
	je	.L40
	.loc 3 323 0
	movq	-24(%rbp), %rax
	movl	820(%rax), %eax
	cmpl	$1, %eax
	je	.L41
	movq	-24(%rbp), %rax
	movq	464(%rax), %rax
	movq	(%rax), %rax
	testq	%rax, %rax
	jne	.L42
.L41:
	.loc 3 324 0
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	VecCUSPCopyFromGPU
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L42
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC4, %ecx
	movl	$__func__.15326, %edx
	movl	$324, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L43
.L42:
	.loc 3 327 0
	movq	-24(%rbp), %rax
	movq	464(%rax), %rax
	movq	(%rax), %rdx
	movq	-32(%rbp), %rax
	movq	%rdx, (%rax)
	jmp	.L44
.L40:
	.loc 3 329 0
	movq	-24(%rbp), %rax
	movq	448(%rax), %rax
	movq	184(%rax), %rcx
	movq	-32(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	*%rcx
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L44
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC4, %ecx
	movl	$__func__.15326, %edx
	movl	$329, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L43
.L44:
	.loc 3 331 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L45
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L45
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L45:
	movl	$0, %eax
.L43:
	.loc 3 332 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE17:
	.size	VecGetArray, .-VecGetArray
	.type	VecRestoreArray, @function
VecRestoreArray:
.LFB18:
	.loc 3 337 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	.loc 3 340 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L59
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L60
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.15408, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC4, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$340, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L58
.L59:
	nop
	jmp	.L58
.L60:
	nop
.L58:
	.loc 3 341 0
	movq	-24(%rbp), %rax
	movl	816(%rax), %eax
	testl	%eax, %eax
	je	.L52
	.loc 3 343 0
	movq	-24(%rbp), %rax
	movl	820(%rax), %eax
	testl	%eax, %eax
	je	.L53
	.loc 3 344 0
	movq	-24(%rbp), %rax
	movl	$2, 820(%rax)
	jmp	.L53
.L52:
	.loc 3 348 0
	movq	-24(%rbp), %rax
	movq	448(%rax), %rax
	movq	208(%rax), %rcx
	movq	-32(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	*%rcx
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L53
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC4, %ecx
	movl	$__func__.15408, %edx
	movl	$348, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L54
.L53:
	.loc 3 350 0
	movq	-24(%rbp), %rax
	movl	164(%rax), %edx
	addl	$1, %edx
	movl	%edx, 164(%rax)
	movl	$0, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L55
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC4, %ecx
	movl	$__func__.15408, %edx
	movl	$350, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L54
.L55:
	.loc 3 351 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L56
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L56
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L56:
	movl	$0, %eax
.L54:
	.loc 3 352 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE18:
	.size	VecRestoreArray, .-VecRestoreArray
	.comm	mtNexts,4,4
	.comm	mtNexti,4,4
	.comm	s_seeds,2496,32
	.section	.rodata
.LC6:
	.string	"baij2.c"
.LC7:
	.string	"Negative overlap specified"
.LC8:
	.string	"index greater than mat-dim"
	.text
.globl MatIncreaseOverlap_SeqBAIJ
	.type	MatIncreaseOverlap_SeqBAIJ, @function
MatIncreaseOverlap_SeqBAIJ:
.LFB62:
	.file 4 "baij2.c"
	.loc 4 10 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$184, %rsp
	movq	%rdi, -152(%rbp)
	movl	%esi, -156(%rbp)
	movq	%rdx, -168(%rbp)
	movl	%ecx, -172(%rbp)
	.loc 4 11 0
	movq	-152(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 18 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L115
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L116
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.24258, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$18, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L114
.L115:
	nop
	jmp	.L114
.L116:
	nop
.L114:
	.loc 4 19 0
	movq	-96(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -64(%rbp)
	.loc 4 20 0
	movq	-96(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -40(%rbp)
	.loc 4 21 0
	movq	-96(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -32(%rbp)
	.loc 4 22 0
	movq	-152(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 24 0
	cmpl	$0, -172(%rbp)
	jns	.L63
	movq	$.LC7, 8(%rsp)
	movl	$0, (%rsp)
	movl	$63, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$24, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L63:
	.loc 4 26 0
	movl	-64(%rbp), %eax
	addl	$15, %eax
	cmpl	$7, %eax
	jbe	.L65
	movq	PetscTrMalloc(%rip), %rbx
	leaq	-136(%rbp), %rdx
	movl	-64(%rbp), %eax
	leal	7(%rax), %ecx
	testl	%eax, %eax
	cmovs	%ecx, %eax
	sarl	$3, %eax
	addl	$1, %eax
	cltq
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$26, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L66
	jmp	.L67
.L65:
	movq	$0, -136(%rbp)
.L67:
	movl	-64(%rbp), %eax
	leal	7(%rax), %edx
	testl	%eax, %eax
	cmovs	%edx, %eax
	sarl	$3, %eax
	addl	$1, %eax
	movslq	%eax, %rdx
	movq	-136(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	testl	%eax, %eax
	je	.L68
.L66:
	movl	$1, %eax
	jmp	.L69
.L68:
	movl	$0, %eax
.L69:
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L70
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$26, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L70:
	.loc 4 27 0
	movl	-64(%rbp), %eax
	addl	$1, %eax
	cltq
	salq	$2, %rax
	testq	%rax, %rax
	je	.L71
	movq	PetscTrMalloc(%rip), %rax
	leaq	-112(%rbp), %rdx
	movl	-64(%rbp), %ecx
	addl	$1, %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,4), %rbx
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$27, %esi
	movq	%rbx, %rdi
	call	*%rax
	jmp	.L72
.L71:
	movq	$0, -112(%rbp)
	movl	$0, %eax
.L72:
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L73
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$27, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L73:
	.loc 4 28 0
	movq	-152(%rbp), %rax
	movq	456(%rax), %rax
	movl	8(%rax), %eax
	addl	$1, %eax
	cltq
	salq	$2, %rax
	testq	%rax, %rax
	je	.L74
	movq	PetscTrMalloc(%rip), %rbx
	leaq	-128(%rbp), %rdx
	movq	-152(%rbp), %rax
	movq	456(%rax), %rax
	movl	8(%rax), %eax
	addl	$1, %eax
	cltq
	salq	$2, %rax
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$28, %esi
	movq	%rax, %rdi
	call	*%rbx
	jmp	.L75
.L74:
	movq	$0, -128(%rbp)
	movl	$0, %eax
.L75:
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L76
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$28, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L76:
	.loc 4 30 0
	movl	$0, -80(%rbp)
	jmp	.L77
.L99:
	.loc 4 32 0
	movl	$0, -60(%rbp)
	.loc 4 33 0
	movl	-64(%rbp), %eax
	leal	7(%rax), %edx
	testl	%eax, %eax
	cmovs	%edx, %eax
	sarl	$3, %eax
	addl	$1, %eax
	movslq	%eax, %rdx
	movq	-136(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L78
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$33, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L78:
	.loc 4 36 0
	movl	-80(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-168(%rbp), %rax
	movq	(%rax), %rax
	leaq	-120(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetIndices
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L79
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$36, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L79:
	.loc 4 37 0
	movl	-80(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-168(%rbp), %rax
	movq	(%rax), %rax
	leaq	-100(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetLocalSize
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L80
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$37, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L80:
	.loc 4 40 0
	movl	$0, -76(%rbp)
	jmp	.L81
.L84:
	.loc 4 41 0
	movq	-120(%rbp), %rax
	movl	-76(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	movl	%eax, %edx
	sarl	$31, %edx
	idivl	-20(%rbp)
	movl	%eax, -52(%rbp)
	.loc 4 42 0
	movl	-52(%rbp), %eax
	cmpl	-64(%rbp), %eax
	jl	.L82
	movq	$.LC8, 8(%rsp)
	movl	$0, (%rsp)
	movl	$63, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$42, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L82:
	.loc 4 43 0
	movl	-52(%rbp), %eax
	leal	7(%rax), %edx
	testl	%eax, %eax
	cmovs	%edx, %eax
	sarl	$3, %eax
	movl	%eax, _BT_idx(%rip)
	movq	-136(%rbp), %rdx
	movl	_BT_idx(%rip), %eax
	cltq
	leaq	(%rdx,%rax), %rax
	movzbl	(%rax), %eax
	movb	%al, _BT_c(%rip)
	movl	-52(%rbp), %eax
	movl	%eax, %edx
	sarl	$31, %edx
	shrl	$29, %edx
	addl	%edx, %eax
	andl	$7, %eax
	subl	%edx, %eax
	movl	$1, %edx
	movl	%edx, %ebx
	movl	%eax, %ecx
	sall	%cl, %ebx
	movl	%ebx, %eax
	movb	%al, _BT_mask(%rip)
	movq	-136(%rbp), %rdx
	movl	_BT_idx(%rip), %eax
	cltq
	addq	%rax, %rdx
	movzbl	_BT_c(%rip), %ecx
	movzbl	_BT_mask(%rip), %eax
	orl	%ecx, %eax
	movb	%al, (%rdx)
	movzbl	_BT_c(%rip), %edx
	movzbl	_BT_mask(%rip), %eax
	andl	%edx, %eax
	testb	%al, %al
	jne	.L83
	movq	-112(%rbp), %rax
	movl	-60(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	leaq	(%rax,%rdx), %rdx
	movl	-52(%rbp), %eax
	movl	%eax, (%rdx)
	addl	$1, -60(%rbp)
.L83:
	.loc 4 40 0
	addl	$1, -76(%rbp)
.L81:
	movl	-100(%rbp), %eax
	cmpl	%eax, -76(%rbp)
	jl	.L84
	.loc 4 45 0
	movl	-80(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-168(%rbp), %rax
	movq	(%rax), %rax
	leaq	-120(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISRestoreIndices
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L85
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$45, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L85:
	.loc 4 46 0
	movl	-80(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-168(%rbp), %rax
	movq	%rax, %rdi
	call	ISDestroy
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L86
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$46, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L86:
	.loc 4 48 0
	movl	$0, -72(%rbp)
	.loc 4 49 0
	movl	$0, -76(%rbp)
	jmp	.L87
.L93:
	.loc 4 50 0
	movl	-60(%rbp), %eax
	movl	%eax, -100(%rbp)
	.loc 4 51 0
	jmp	.L88
.L92:
	.loc 4 52 0
	movq	-112(%rbp), %rax
	movl	-72(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	movl	%eax, -84(%rbp)
	.loc 4 53 0
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-40(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -48(%rbp)
	.loc 4 54 0
	movl	-84(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-40(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -44(%rbp)
	.loc 4 55 0
	movl	-48(%rbp), %eax
	movl	%eax, -68(%rbp)
	jmp	.L89
.L91:
	.loc 4 56 0
	movl	-68(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -56(%rbp)
	.loc 4 57 0
	movl	-56(%rbp), %eax
	leal	7(%rax), %edx
	testl	%eax, %eax
	cmovs	%edx, %eax
	sarl	$3, %eax
	movl	%eax, _BT_idx(%rip)
	movq	-136(%rbp), %rdx
	movl	_BT_idx(%rip), %eax
	cltq
	leaq	(%rdx,%rax), %rax
	movzbl	(%rax), %eax
	movb	%al, _BT_c(%rip)
	movl	-56(%rbp), %eax
	movl	%eax, %edx
	sarl	$31, %edx
	shrl	$29, %edx
	addl	%edx, %eax
	andl	$7, %eax
	subl	%edx, %eax
	movl	$1, %edx
	movl	%edx, %ebx
	movl	%eax, %ecx
	sall	%cl, %ebx
	movl	%ebx, %eax
	movb	%al, _BT_mask(%rip)
	movq	-136(%rbp), %rdx
	movl	_BT_idx(%rip), %eax
	cltq
	addq	%rax, %rdx
	movzbl	_BT_c(%rip), %ecx
	movzbl	_BT_mask(%rip), %eax
	orl	%ecx, %eax
	movb	%al, (%rdx)
	movzbl	_BT_c(%rip), %edx
	movzbl	_BT_mask(%rip), %eax
	andl	%edx, %eax
	testb	%al, %al
	jne	.L90
	movq	-112(%rbp), %rax
	movl	-60(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	leaq	(%rax,%rdx), %rdx
	movl	-56(%rbp), %eax
	movl	%eax, (%rdx)
	addl	$1, -60(%rbp)
.L90:
	.loc 4 55 0
	addl	$1, -68(%rbp)
.L89:
	movl	-68(%rbp), %eax
	cmpl	-44(%rbp), %eax
	jl	.L91
	.loc 4 51 0
	addl	$1, -72(%rbp)
.L88:
	movl	-100(%rbp), %eax
	cmpl	%eax, -72(%rbp)
	jl	.L92
	.loc 4 49 0
	addl	$1, -76(%rbp)
.L87:
	movl	-76(%rbp), %eax
	cmpl	-172(%rbp), %eax
	jl	.L93
	.loc 4 62 0
	movl	$0, -76(%rbp)
	jmp	.L94
.L97:
	.loc 4 63 0
	movl	$0, -72(%rbp)
	jmp	.L95
.L96:
	.loc 4 64 0
	movq	-128(%rbp), %rdx
	movl	-76(%rbp), %eax
	imull	-20(%rbp), %eax
	addl	-72(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	%rax, %rdx
	movq	-112(%rbp), %rax
	movl	-76(%rbp), %ecx
	movslq	%ecx, %rcx
	salq	$2, %rcx
	addq	%rcx, %rax
	movl	(%rax), %eax
	imull	-20(%rbp), %eax
	addl	-72(%rbp), %eax
	movl	%eax, (%rdx)
	.loc 4 63 0
	addl	$1, -72(%rbp)
.L95:
	movl	-72(%rbp), %eax
	cmpl	-20(%rbp), %eax
	jl	.L96
	.loc 4 62 0
	addl	$1, -76(%rbp)
.L94:
	movl	-76(%rbp), %eax
	cmpl	-60(%rbp), %eax
	jl	.L97
	.loc 4 66 0
	movl	-80(%rbp), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	addq	-168(%rbp), %rcx
	movq	-128(%rbp), %rdx
	movl	-60(%rbp), %eax
	imull	-20(%rbp), %eax
	movq	%rcx, %r8
	movl	$0, %ecx
	movl	%eax, %esi
	movl	$1140850689, %edi
	call	ISCreateGeneral
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L98
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$66, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L98:
	.loc 4 30 0
	addl	$1, -80(%rbp)
.L77:
	movl	-80(%rbp), %eax
	cmpl	-156(%rbp), %eax
	jl	.L99
	.loc 4 68 0
	movq	-136(%rbp), %rax
	testq	%rax, %rax
	je	.L100
	movq	PetscTrFree(%rip), %rbx
	movq	-136(%rbp), %rax
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$68, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L101
	movq	$0, -136(%rbp)
	jmp	.L100
.L101:
	movl	$1, %eax
	jmp	.L102
.L100:
	movl	$0, %eax
.L102:
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L103
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$68, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L103:
	.loc 4 69 0
	movq	-112(%rbp), %rax
	testq	%rax, %rax
	je	.L104
	movq	PetscTrFree(%rip), %rbx
	movq	-112(%rbp), %rax
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$69, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L105
	movq	$0, -112(%rbp)
	jmp	.L104
.L105:
	movl	$1, %eax
	jmp	.L106
.L104:
	movl	$0, %eax
.L106:
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L107
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$69, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L107:
	.loc 4 70 0
	movq	-128(%rbp), %rax
	testq	%rax, %rax
	je	.L108
	movq	PetscTrFree(%rip), %rbx
	movq	-128(%rbp), %rax
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$70, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L109
	movq	$0, -128(%rbp)
	jmp	.L108
.L109:
	movl	$1, %eax
	jmp	.L110
.L108:
	movl	$0, %eax
.L110:
	movl	%eax, -88(%rbp)
	cmpl	$0, -88(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L111
	movl	-88(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24258, %edx
	movl	$70, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L64
.L111:
	.loc 4 71 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L112
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L112
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L112:
	movl	$0, %eax
.L64:
	.loc 4 72 0
	addq	$184, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE62:
	.size	MatIncreaseOverlap_SeqBAIJ, .-MatIncreaseOverlap_SeqBAIJ
	.section	.rodata
.LC9:
	.string	"IS is not sorted"
.LC10:
	.string	"Submatrix wrong size"
	.align 8
.LC11:
	.string	"Cannot reuse matrix. wrong no of nonzeros"
	.text
.globl MatGetSubMatrix_SeqBAIJ_Private
	.type	MatGetSubMatrix_SeqBAIJ_Private, @function
MatGetSubMatrix_SeqBAIJ_Private:
.LFB63:
	.loc 4 77 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$248, %rsp
	movq	%rdi, -200(%rbp)
	movq	%rsi, -208(%rbp)
	movq	%rdx, -216(%rbp)
	movl	%ecx, -220(%rbp)
	movq	%r8, -232(%rbp)
	.loc 4 78 0
	movq	-200(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -128(%rbp)
	.loc 4 80 0
	movq	-128(%rbp), %rax
	movl	232(%rax), %eax
	movl	%eax, -92(%rbp)
	.loc 4 83 0
	movq	-200(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -48(%rbp)
	movq	-128(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -44(%rbp)
	.loc 4 84 0
	movq	-128(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -40(%rbp)
	movq	-128(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -32(%rbp)
	.loc 4 89 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L171
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L172
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.24627, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$89, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L170
.L171:
	nop
	jmp	.L170
.L172:
	nop
.L170:
	.loc 4 90 0
	leaq	-184(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISSorted
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L119
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$90, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L119:
	.loc 4 91 0
	movl	-184(%rbp), %eax
	testl	%eax, %eax
	jne	.L121
	movq	$.LC9, 8(%rsp)
	movl	$0, (%rsp)
	movl	$73, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$91, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L121:
	.loc 4 93 0
	leaq	-152(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetIndices
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L122
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$93, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L122:
	.loc 4 94 0
	leaq	-160(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetIndices
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L123
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$94, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L123:
	.loc 4 95 0
	leaq	-164(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetLocalSize
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L124
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$95, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L124:
	.loc 4 96 0
	leaq	-168(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetLocalSize
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L125
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$96, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L125:
	.loc 4 98 0
	movl	-92(%rbp), %eax
	addl	$1, %eax
	cltq
	salq	$2, %rax
	testq	%rax, %rax
	je	.L126
	movq	PetscTrMalloc(%rip), %rax
	leaq	-136(%rbp), %rdx
	movl	-92(%rbp), %ecx
	addl	$1, %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,4), %rbx
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$98, %esi
	movq	%rbx, %rdi
	call	*%rax
	jmp	.L127
.L126:
	movq	$0, -136(%rbp)
	movl	$0, %eax
.L127:
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L128
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$98, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L128:
	.loc 4 99 0
	movq	-136(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 100 0
	movl	-164(%rbp), %eax
	addl	$1, %eax
	cltq
	salq	$2, %rax
	testq	%rax, %rax
	je	.L129
	movq	PetscTrMalloc(%rip), %rax
	leaq	-144(%rbp), %rdx
	movl	-164(%rbp), %ecx
	addl	$1, %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,4), %rbx
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$100, %esi
	movq	%rbx, %rdi
	call	*%rax
	jmp	.L130
.L129:
	movq	$0, -144(%rbp)
	movl	$0, %eax
.L130:
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L131
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$100, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L131:
	.loc 4 101 0
	movl	-92(%rbp), %eax
	cltq
	leaq	0(,%rax,4), %rdx
	movq	-136(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L132
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$101, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L132:
	.loc 4 102 0
	movl	$0, -108(%rbp)
	jmp	.L133
.L134:
	movq	-136(%rbp), %rax
	movq	-160(%rbp), %rdx
	movl	-108(%rbp), %ecx
	movslq	%ecx, %rcx
	salq	$2, %rcx
	addq	%rcx, %rdx
	movl	(%rdx), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	-108(%rbp), %edx
	addl	$1, %edx
	movl	%edx, (%rax)
	addl	$1, -108(%rbp)
.L133:
	movl	-168(%rbp), %eax
	cmpl	%eax, -108(%rbp)
	jl	.L134
	.loc 4 104 0
	movl	$0, -108(%rbp)
	jmp	.L135
.L139:
	.loc 4 105 0
	movq	-152(%rbp), %rax
	movl	-108(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -100(%rbp)
	.loc 4 106 0
	movq	-128(%rbp), %rax
	movq	32(%rax), %rax
	movq	-152(%rbp), %rdx
	movl	-108(%rbp), %ecx
	movslq	%ecx, %rcx
	salq	$2, %rcx
	addq	%rcx, %rdx
	movl	(%rdx), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	addl	-100(%rbp), %eax
	movl	%eax, -96(%rbp)
	.loc 4 107 0
	movq	-144(%rbp), %rax
	movl	-108(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	$0, (%rax)
	.loc 4 108 0
	movl	-100(%rbp), %eax
	movl	%eax, -104(%rbp)
	jmp	.L136
.L138:
	.loc 4 109 0
	movl	-104(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-40(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	testl	%eax, %eax
	je	.L137
	.loc 4 110 0
	movq	-144(%rbp), %rax
	movl	-108(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %edx
	addl	$1, %edx
	movl	%edx, (%rax)
.L137:
	.loc 4 108 0
	addl	$1, -104(%rbp)
.L136:
	movl	-104(%rbp), %eax
	cmpl	-96(%rbp), %eax
	jl	.L138
	.loc 4 104 0
	addl	$1, -108(%rbp)
.L135:
	movl	-164(%rbp), %eax
	cmpl	%eax, -108(%rbp)
	jl	.L139
	.loc 4 115 0
	cmpl	$1, -220(%rbp)
	jne	.L140
	.loc 4 116 0
	movq	-232(%rbp), %rax
	movq	(%rax), %rax
	movq	472(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 118 0
	movq	-120(%rbp), %rax
	movl	228(%rax), %edx
	movl	-164(%rbp), %eax
	cmpl	%eax, %edx
	jne	.L141
	movq	-120(%rbp), %rax
	movl	232(%rax), %edx
	movl	-168(%rbp), %eax
	cmpl	%eax, %edx
	jne	.L141
	movq	-232(%rbp), %rax
	movq	(%rax), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	cmpl	-48(%rbp), %eax
	je	.L142
.L141:
	movq	$.LC10, 8(%rsp)
	movl	$0, (%rsp)
	movl	$60, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$118, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L142:
	.loc 4 119 0
	movq	-120(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	leaq	0(,%rax,4), %rsi
	movq	-144(%rbp), %rbx
	movq	-120(%rbp), %rax
	movq	32(%rax), %rax
	leaq	-180(%rbp), %rdx
	movq	%rdx, %rcx
	movq	%rsi, %rdx
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	PetscMemcmp
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L143
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$119, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L143:
	.loc 4 120 0
	movl	-180(%rbp), %eax
	testl	%eax, %eax
	jne	.L144
	movq	$.LC11, 8(%rsp)
	movl	$0, (%rsp)
	movl	$60, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$120, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L144:
	.loc 4 121 0
	movq	-120(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	leaq	0(,%rax,4), %rdx
	movq	-120(%rbp), %rax
	movq	32(%rax), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L145
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$121, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L145:
	.loc 4 122 0
	movq	-232(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	jmp	.L146
.L140:
	.loc 4 124 0
	movq	-200(%rbp), %rax
	movl	16(%rax), %eax
	leaq	-176(%rbp), %rdx
	movq	%rdx, %rsi
	movl	%eax, %edi
	call	MatCreate
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L147
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$124, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L147:
	.loc 4 125 0
	movl	-168(%rbp), %eax
	movl	%eax, %edx
	imull	-48(%rbp), %edx
	movl	-164(%rbp), %eax
	movl	%eax, %ebx
	imull	-48(%rbp), %ebx
	movq	-176(%rbp), %rax
	movl	$-1, %r8d
	movl	$-1, %ecx
	movl	%ebx, %esi
	movq	%rax, %rdi
	call	MatSetSizes
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L148
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$125, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L148:
	.loc 4 126 0
	movq	-200(%rbp), %rax
	movq	104(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	MatSetType
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L149
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$126, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L149:
	.loc 4 127 0
	movq	-144(%rbp), %rdx
	movq	-176(%rbp), %rax
	movl	-48(%rbp), %ebx
	movq	%rdx, %rcx
	movl	$0, %edx
	movl	%ebx, %esi
	movq	%rax, %rdi
	call	MatSeqBAIJSetPreallocation_SeqBAIJ
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L146
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$127, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L146:
	.loc 4 129 0
	movq	-176(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 130 0
	movl	$0, -108(%rbp)
	jmp	.L150
.L155:
	.loc 4 131 0
	movq	-152(%rbp), %rax
	movl	-108(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	movl	%eax, -88(%rbp)
	.loc 4 132 0
	movl	-88(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -100(%rbp)
	.loc 4 133 0
	movq	-128(%rbp), %rax
	movq	32(%rax), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	addl	-100(%rbp), %eax
	movl	%eax, -96(%rbp)
	.loc 4 134 0
	movq	-120(%rbp), %rax
	movq	136(%rax), %rax
	movl	-108(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	movl	%eax, -84(%rbp)
	.loc 4 135 0
	movq	-120(%rbp), %rax
	movq	144(%rax), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -80(%rbp)
	.loc 4 136 0
	movq	-120(%rbp), %rax
	movq	168(%rax), %rdx
	movl	-84(%rbp), %eax
	imull	-44(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -24(%rbp)
	.loc 4 137 0
	movq	-120(%rbp), %rax
	movq	32(%rax), %rax
	movl	-108(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -64(%rbp)
	.loc 4 138 0
	movl	-100(%rbp), %eax
	movl	%eax, -104(%rbp)
	jmp	.L151
.L154:
	.loc 4 139 0
	movq	-128(%rbp), %rax
	movq	144(%rax), %rax
	movl	-104(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -68(%rbp)
	cmpl	$0, -68(%rbp)
	je	.L152
	.loc 4 140 0
	movl	-68(%rbp), %eax
	leal	-1(%rax), %edx
	movq	-80(%rbp), %rax
	movl	%edx, (%rax)
	addq	$4, -80(%rbp)
	.loc 4 141 0
	movl	-44(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-128(%rbp), %rax
	movq	168(%rax), %rcx
	movl	-104(%rbp), %eax
	imull	-44(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, %rcx
	movq	-24(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L153
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$141, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L153:
	.loc 4 142 0
	movl	-44(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -24(%rbp)
	.loc 4 143 0
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	leal	1(%rax), %edx
	movq	-64(%rbp), %rax
	movl	%edx, (%rax)
.L152:
	.loc 4 138 0
	addl	$1, -104(%rbp)
.L151:
	movl	-104(%rbp), %eax
	cmpl	-96(%rbp), %eax
	jl	.L154
	.loc 4 130 0
	addl	$1, -108(%rbp)
.L150:
	movl	-164(%rbp), %eax
	cmpl	%eax, -108(%rbp)
	jl	.L155
	.loc 4 149 0
	leaq	-160(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISRestoreIndices
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L156
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$149, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L156:
	.loc 4 150 0
	movq	-136(%rbp), %rax
	testq	%rax, %rax
	je	.L157
	movq	PetscTrFree(%rip), %rbx
	movq	-136(%rbp), %rax
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$150, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L158
	movq	$0, -136(%rbp)
	jmp	.L157
.L158:
	movl	$1, %eax
	jmp	.L159
.L157:
	movl	$0, %eax
.L159:
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L160
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$150, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L160:
	.loc 4 151 0
	movq	-144(%rbp), %rax
	testq	%rax, %rax
	je	.L161
	movq	PetscTrFree(%rip), %rbx
	movq	-144(%rbp), %rax
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$151, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L162
	movq	$0, -144(%rbp)
	jmp	.L161
.L162:
	movl	$1, %eax
	jmp	.L163
.L161:
	movl	$0, %eax
.L163:
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L164
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$151, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L164:
	.loc 4 152 0
	movq	-176(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	MatAssemblyBegin
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L165
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$152, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L165:
	.loc 4 153 0
	movq	-176(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	MatAssemblyEnd
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L166
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$153, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L166:
	.loc 4 155 0
	leaq	-152(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISRestoreIndices
	movl	%eax, -112(%rbp)
	cmpl	$0, -112(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L167
	movl	-112(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.24627, %edx
	movl	$155, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L120
.L167:
	.loc 4 156 0
	movq	-176(%rbp), %rdx
	movq	-232(%rbp), %rax
	movq	%rdx, (%rax)
	.loc 4 157 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L168
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L168
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L168:
	movl	$0, %eax
.L120:
	.loc 4 158 0
	addq	$248, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE63:
	.size	MatGetSubMatrix_SeqBAIJ_Private, .-MatGetSubMatrix_SeqBAIJ_Private
	.section	.rodata
	.align 8
.LC12:
	.string	"Index set does not match blocks"
.LC13:
	.string	"Internal error in PETSc"
	.text
.globl MatGetSubMatrix_SeqBAIJ
	.type	MatGetSubMatrix_SeqBAIJ, @function
MatGetSubMatrix_SeqBAIJ:
.LFB64:
	.loc 4 163 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$152, %rsp
	movq	%rdi, -104(%rbp)
	movq	%rsi, -112(%rbp)
	movq	%rdx, -120(%rbp)
	movl	%ecx, -124(%rbp)
	movq	%r8, -136(%rbp)
	.loc 4 164 0
	movq	-104(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -40(%rbp)
	.loc 4 167 0
	movq	-104(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -24(%rbp)
	.loc 4 170 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L213
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L214
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.25037, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$170, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L212
.L213:
	nop
	jmp	.L212
.L214:
	nop
.L212:
	.loc 4 171 0
	leaq	-88(%rbp), %rdx
	movq	-112(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetIndices
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L175
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$171, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L175:
	.loc 4 172 0
	leaq	-96(%rbp), %rdx
	movq	-120(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetIndices
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L177
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$172, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L177:
	.loc 4 173 0
	leaq	-76(%rbp), %rdx
	movq	-112(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetLocalSize
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L178
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$173, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L178:
	.loc 4 174 0
	leaq	-80(%rbp), %rdx
	movq	-120(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISGetLocalSize
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L179
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$174, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L179:
	.loc 4 178 0
	movq	$0, -72(%rbp)
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	salq	$3, %rax
	cmpq	$-15, %rax
	je	.L180
	movq	PetscTrMalloc(%rip), %rbx
	leaq	-64(%rbp), %rdx
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	salq	$3, %rax
	addq	$15, %rax
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$178, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	setne	%al
	jmp	.L181
.L180:
	movq	$0, -64(%rbp)
	movl	$0, %eax
.L181:
	testb	%al, %al
	jne	.L182
	movq	-64(%rbp), %rdx
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	addq	$15, %rax
	andq	$-16, %rax
	movq	%rax, -72(%rbp)
	movl	$0, %eax
	jmp	.L183
.L182:
	movl	$1, %eax
.L183:
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L184
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$178, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L184:
	.loc 4 179 0
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	leaq	0(,%rax,4), %rdx
	movq	-64(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L185
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$179, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L185:
	.loc 4 180 0
	movl	$0, -28(%rbp)
	jmp	.L186
.L187:
	movq	-64(%rbp), %rcx
	movq	-88(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	movl	%eax, %edx
	sarl	$31, %edx
	idivl	-24(%rbp)
	cltq
	salq	$2, %rax
	leaq	(%rcx,%rax), %rax
	movl	(%rax), %edx
	addl	$1, %edx
	movl	%edx, (%rax)
	addl	$1, -28(%rbp)
.L186:
	movl	-76(%rbp), %eax
	cmpl	%eax, -28(%rbp)
	jl	.L187
	.loc 4 181 0
	movl	$0, -20(%rbp)
	.loc 4 182 0
	movl	$0, -28(%rbp)
	jmp	.L188
.L191:
	.loc 4 183 0
	movq	-64(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	testl	%eax, %eax
	je	.L189
	movq	-64(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	cmpl	-24(%rbp), %eax
	je	.L189
	movq	$.LC12, 8(%rsp)
	movl	$0, (%rsp)
	movl	$60, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$183, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L189:
	.loc 4 184 0
	movq	-64(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	cmpl	-24(%rbp), %eax
	jne	.L190
	movq	-72(%rbp), %rax
	movl	-20(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	leaq	(%rax,%rdx), %rdx
	movl	-28(%rbp), %eax
	movl	%eax, (%rdx)
	addl	$1, -20(%rbp)
.L190:
	.loc 4 182 0
	addl	$1, -28(%rbp)
.L188:
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	cmpl	-28(%rbp), %eax
	jg	.L191
	.loc 4 186 0
	movq	-72(%rbp), %rdx
	leaq	-48(%rbp), %rcx
	movl	-20(%rbp), %eax
	movq	%rcx, %r8
	movl	$0, %ecx
	movl	%eax, %esi
	movl	$1140850689, %edi
	call	ISCreateGeneral
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L192
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$186, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L192:
	.loc 4 188 0
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	leaq	0(,%rax,4), %rdx
	movq	-64(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L193
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$188, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L193:
	.loc 4 189 0
	movl	$0, -28(%rbp)
	jmp	.L194
.L195:
	movq	-64(%rbp), %rcx
	movq	-96(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	movl	%eax, %edx
	sarl	$31, %edx
	idivl	-24(%rbp)
	cltq
	salq	$2, %rax
	leaq	(%rcx,%rax), %rax
	movl	(%rax), %edx
	addl	$1, %edx
	movl	%edx, (%rax)
	addl	$1, -28(%rbp)
.L194:
	movl	-80(%rbp), %eax
	cmpl	%eax, -28(%rbp)
	jl	.L195
	.loc 4 190 0
	movl	$0, -20(%rbp)
	.loc 4 191 0
	movl	$0, -28(%rbp)
	jmp	.L196
.L199:
	.loc 4 192 0
	movq	-64(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	testl	%eax, %eax
	je	.L197
	movq	-64(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	cmpl	-24(%rbp), %eax
	je	.L197
	movq	$.LC13, 8(%rsp)
	movl	$0, (%rsp)
	movl	$77, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$192, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L197:
	.loc 4 193 0
	movq	-64(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %eax
	cmpl	-24(%rbp), %eax
	jne	.L198
	movq	-72(%rbp), %rax
	movl	-20(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	leaq	(%rax,%rdx), %rdx
	movl	-28(%rbp), %eax
	movl	%eax, (%rdx)
	addl	$1, -20(%rbp)
.L198:
	.loc 4 191 0
	addl	$1, -28(%rbp)
.L196:
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	cmpl	-28(%rbp), %eax
	jg	.L199
	.loc 4 195 0
	movq	-72(%rbp), %rdx
	leaq	-56(%rbp), %rcx
	movl	-20(%rbp), %eax
	movq	%rcx, %r8
	movl	$0, %ecx
	movl	%eax, %esi
	movl	$1140850689, %edi
	call	ISCreateGeneral
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L200
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$195, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L200:
	.loc 4 196 0
	leaq	-88(%rbp), %rdx
	movq	-112(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISRestoreIndices
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L201
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$196, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L201:
	.loc 4 197 0
	leaq	-96(%rbp), %rdx
	movq	-120(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	ISRestoreIndices
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L202
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$197, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L202:
	.loc 4 198 0
	movq	$0, -72(%rbp)
	movq	-64(%rbp), %rax
	testq	%rax, %rax
	je	.L203
	movq	PetscTrFree(%rip), %rbx
	movq	-64(%rbp), %rax
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$198, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L204
	movq	$0, -64(%rbp)
	jmp	.L203
.L204:
	movl	$1, %eax
	jmp	.L205
.L203:
	movl	$0, %eax
.L205:
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L206
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$198, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L206:
	.loc 4 200 0
	movq	-56(%rbp), %rdx
	movq	-48(%rbp), %rbx
	movq	-136(%rbp), %rsi
	movl	-124(%rbp), %ecx
	movq	-104(%rbp), %rax
	movq	%rsi, %r8
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	MatGetSubMatrix_SeqBAIJ_Private
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L207
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$200, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L207:
	.loc 4 201 0
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	ISDestroy
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L208
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$201, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L208:
	.loc 4 202 0
	leaq	-56(%rbp), %rax
	movq	%rax, %rdi
	call	ISDestroy
	movl	%eax, -32(%rbp)
	cmpl	$0, -32(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L209
	movl	-32(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25037, %edx
	movl	$202, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L176
.L209:
	.loc 4 203 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L210
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L210
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L210:
	movl	$0, %eax
.L176:
	.loc 4 204 0
	addq	$152, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE64:
	.size	MatGetSubMatrix_SeqBAIJ, .-MatGetSubMatrix_SeqBAIJ
.globl MatGetSubMatrices_SeqBAIJ
	.type	MatGetSubMatrices_SeqBAIJ, @function
MatGetSubMatrices_SeqBAIJ:
.LFB65:
	.loc 4 209 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$88, %rsp
	movq	%rdi, -40(%rbp)
	movl	%esi, -44(%rbp)
	movq	%rdx, -56(%rbp)
	movq	%rcx, -64(%rbp)
	movl	%r8d, -68(%rbp)
	movq	%r9, -80(%rbp)
	.loc 4 213 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L227
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L228
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.25328, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$213, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L226
.L227:
	nop
	jmp	.L226
.L228:
	nop
.L226:
	.loc 4 214 0
	cmpl	$0, -68(%rbp)
	jne	.L217
	.loc 4 215 0
	movl	-44(%rbp), %eax
	addl	$1, %eax
	cltq
	salq	$3, %rax
	testq	%rax, %rax
	je	.L218
	movq	PetscTrMalloc(%rip), %rax
	movq	-80(%rbp), %rdx
	movl	-44(%rbp), %ecx
	addl	$1, %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,8), %rbx
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25328, %edx
	movl	$215, %esi
	movq	%rbx, %rdi
	call	*%rax
	jmp	.L219
.L218:
	movq	-80(%rbp), %rax
	movq	$0, (%rax)
	movl	$0, %eax
.L219:
	movl	%eax, -24(%rbp)
	cmpl	$0, -24(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L217
	movl	-24(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25328, %edx
	movl	$215, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L220
.L217:
	.loc 4 218 0
	movl	$0, -20(%rbp)
	jmp	.L221
.L223:
	.loc 4 219 0
	movq	-80(%rbp), %rax
	movq	(%rax), %rax
	movl	-20(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rsi
	movl	-20(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-64(%rbp), %rax
	movq	(%rax), %rdx
	movl	-20(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-56(%rbp), %rax
	movq	(%rax), %rbx
	movl	-68(%rbp), %ecx
	movq	-40(%rbp), %rax
	movq	%rsi, %r8
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	MatGetSubMatrix_SeqBAIJ
	movl	%eax, -24(%rbp)
	cmpl	$0, -24(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L222
	movl	-24(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25328, %edx
	movl	$219, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L220
.L222:
	.loc 4 218 0
	addl	$1, -20(%rbp)
.L221:
	movl	-20(%rbp), %eax
	cmpl	-44(%rbp), %eax
	jl	.L223
	.loc 4 221 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L224
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L224
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L224:
	movl	$0, %eax
.L220:
	.loc 4 222 0
	addq	$88, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE65:
	.size	MatGetSubMatrices_SeqBAIJ, .-MatGetSubMatrices_SeqBAIJ
.globl MatMult_SeqBAIJ_1
	.type	MatMult_SeqBAIJ_1, @function
MatMult_SeqBAIJ_1:
.LFB66:
	.loc 4 232 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$192, %rsp
	movq	%rdi, -152(%rbp)
	movq	%rsi, -160(%rbp)
	movq	%rdx, -168(%rbp)
	.loc 4 233 0
	movq	-152(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 238 0
	movl	$0, -76(%rbp)
	.loc 4 239 0
	movq	$0, -56(%rbp)
	.loc 4 240 0
	movq	-120(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -44(%rbp)
	.loc 4 242 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L252
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L253
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.25437, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$242, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L251
.L252:
	nop
	jmp	.L251
.L253:
	nop
.L251:
	.loc 4 243 0
	leaq	-136(%rbp), %rdx
	movq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -92(%rbp)
	cmpl	$0, -92(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L231
	movl	-92(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25437, %edx
	movl	$243, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L232
.L231:
	.loc 4 244 0
	leaq	-128(%rbp), %rdx
	movq	-168(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -92(%rbp)
	cmpl	$0, -92(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L233
	movl	-92(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25437, %edx
	movl	$244, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L232
.L233:
	.loc 4 246 0
	cmpl	$0, -44(%rbp)
	je	.L234
	.loc 4 247 0
	movq	-120(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -88(%rbp)
	.loc 4 248 0
	movq	-120(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 249 0
	movq	-120(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -56(%rbp)
	.loc 4 250 0
	movl	-88(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-128(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -92(%rbp)
	cmpl	$0, -92(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L235
	movl	-92(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25437, %edx
	movl	$250, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L232
.L234:
	.loc 4 252 0
	movq	-120(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -88(%rbp)
	.loc 4 253 0
	movq	-120(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
.L235:
	.loc 4 256 0
	movl	$0, -84(%rbp)
	jmp	.L236
.L245:
	.loc 4 257 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -80(%rbp)
	.loc 4 258 0
	movq	-120(%rbp), %rax
	movq	168(%rax), %rdx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 259 0
	movq	-120(%rbp), %rax
	movq	144(%rax), %rdx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 260 0
	addq	$4, -64(%rbp)
.LBB2:
	.loc 4 261 0
	movq	-72(%rbp), %rax
	movl	-80(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -40(%rbp)
	movl	-80(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -32(%rbp)
	jmp	.L237
.L238:
	addq	$64, -40(%rbp)
.L237:
	movq	-40(%rbp), %rax
	cmpq	-32(%rbp), %rax
	jb	.L238
.LBE2:
.LBB3:
	.loc 4 262 0
	movq	-104(%rbp), %rax
	movl	-80(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	%rax, -24(%rbp)
	movl	-80(%rbp), %eax
	cltq
	salq	$4, %rax
	addq	-104(%rbp), %rax
	movq	%rax, -16(%rbp)
	jmp	.L239
.L240:
	addq	$64, -24(%rbp)
.L239:
	movq	-24(%rbp), %rax
	cmpq	-16(%rbp), %rax
	jb	.L240
.LBE3:
	.loc 4 263 0
	movl	$0, %eax
	movq	%rax, -112(%rbp)
.LBB4:
	.loc 4 264 0
	movl	$0, -4(%rbp)
	jmp	.L241
.L242:
	movl	-4(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-104(%rbp), %rax
	vmovsd	(%rax), %xmm1
	movq	-136(%rbp), %rdx
	movl	-4(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	-112(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -112(%rbp)
	addl	$1, -4(%rbp)
.L241:
	movl	-4(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L242
.LBE4:
	.loc 4 265 0
	cmpl	$0, -44(%rbp)
	je	.L243
	.loc 4 266 0
	movq	-128(%rbp), %rdx
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-112(%rbp), %rax
	movq	%rax, (%rdx)
	jmp	.L244
.L243:
	.loc 4 268 0
	cmpl	$0, -80(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -76(%rbp)
	.loc 4 269 0
	movq	-128(%rbp), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rdx
	movq	-112(%rbp), %rax
	movq	%rax, (%rdx)
.L244:
	.loc 4 256 0
	addl	$1, -84(%rbp)
.L236:
	movl	-84(%rbp), %eax
	cmpl	-88(%rbp), %eax
	jl	.L245
	.loc 4 272 0
	leaq	-136(%rbp), %rdx
	movq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -92(%rbp)
	cmpl	$0, -92(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L246
	movl	-92(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25437, %edx
	movl	$272, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L232
.L246:
	.loc 4 273 0
	leaq	-128(%rbp), %rdx
	movq	-168(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -92(%rbp)
	cmpl	$0, -92(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L247
	movl	-92(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25437, %edx
	movl	$273, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L232
.L247:
	.loc 4 274 0
	movq	-120(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	vcvtsi2sd	-76(%rbp), %xmm1, %xmm1
	vsubsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -92(%rbp)
	cmpl	$0, -92(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L248
	movl	-92(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25437, %edx
	movl	$274, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L232
.L248:
	.loc 4 275 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L249
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L249
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L249:
	movl	$0, %eax
.L232:
	.loc 4 276 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE66:
	.size	MatMult_SeqBAIJ_1, .-MatMult_SeqBAIJ_1
.globl MatMult_SeqBAIJ_2
	.type	MatMult_SeqBAIJ_2, @function
MatMult_SeqBAIJ_2:
.LFB67:
	.loc 4 281 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$224, %rsp
	movq	%rdi, -184(%rbp)
	movq	%rsi, -192(%rbp)
	movq	%rdx, -200(%rbp)
	.loc 4 282 0
	movq	-184(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -152(%rbp)
	.loc 4 283 0
	movq	$0, -144(%rbp)
	.loc 4 288 0
	movq	$0, -48(%rbp)
	movl	$0, -40(%rbp)
	.loc 4 289 0
	movq	-152(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 291 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L277
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L278
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.25630, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$291, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L276
.L277:
	nop
	jmp	.L276
.L278:
	nop
.L276:
	.loc 4 292 0
	leaq	-168(%rbp), %rdx
	movq	-192(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L256
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25630, %edx
	movl	$292, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L257
.L256:
	.loc 4 293 0
	leaq	-160(%rbp), %rdx
	movq	-200(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L258
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25630, %edx
	movl	$293, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L257
.L258:
	.loc 4 295 0
	movq	-152(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 296 0
	movq	-152(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 297 0
	cmpl	$0, -36(%rbp)
	je	.L259
	.loc 4 298 0
	movq	-152(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 299 0
	movq	-152(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 300 0
	movq	-152(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L260
.L259:
	.loc 4 302 0
	movq	-152(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 303 0
	movq	-152(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 304 0
	movq	-160(%rbp), %rax
	movq	%rax, -144(%rbp)
.L260:
	.loc 4 307 0
	movl	$0, -76(%rbp)
	jmp	.L261
.L270:
	.loc 4 308 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 309 0
	movl	$0, %eax
	movq	%rax, -136(%rbp)
	movl	$0, %eax
	movq	%rax, -128(%rbp)
	.loc 4 310 0
	cmpl	$0, -52(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -40(%rbp)
.LBB5:
	.loc 4 311 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L262
.L263:
	addq	$64, -32(%rbp)
.L262:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L263
.LBE5:
.LBB6:
	.loc 4 312 0
	movq	-96(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$5, %rdx
	addq	%rdx, %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$6, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L264
.L265:
	addq	$64, -16(%rbp)
.L264:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L265
.LBE6:
	.loc 4 313 0
	movl	$0, -56(%rbp)
	jmp	.L266
.L267:
	.loc 4 314 0
	movq	-168(%rbp), %rdx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$4, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -120(%rbp)
	addq	$4, -72(%rbp)
	movq	-120(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-120(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 315 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 316 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-128(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -128(%rbp)
	.loc 4 317 0
	addq	$32, -96(%rbp)
	.loc 4 313 0
	addl	$1, -56(%rbp)
.L266:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L267
	.loc 4 319 0
	cmpl	$0, -36(%rbp)
	je	.L268
	movq	-160(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$4, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -144(%rbp)
.L268:
	.loc 4 320 0
	movq	-144(%rbp), %rax
	movq	-136(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-144(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-128(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 321 0
	cmpl	$0, -36(%rbp)
	jne	.L269
	addq	$16, -144(%rbp)
.L269:
	.loc 4 307 0
	addl	$1, -76(%rbp)
.L261:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L270
	.loc 4 323 0
	leaq	-168(%rbp), %rdx
	movq	-192(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L271
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25630, %edx
	movl	$323, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L257
.L271:
	.loc 4 324 0
	leaq	-160(%rbp), %rdx
	movq	-200(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L272
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25630, %edx
	movl	$324, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L257
.L272:
	.loc 4 325 0
	movq	-152(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC14(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-40(%rbp), %xmm0, %xmm0
	vmovsd	.LC15(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L273
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25630, %edx
	movl	$325, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L257
.L273:
	.loc 4 326 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L274
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L274
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L274:
	movl	$0, %eax
.L257:
	.loc 4 327 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE67:
	.size	MatMult_SeqBAIJ_2, .-MatMult_SeqBAIJ_2
.globl MatMult_SeqBAIJ_3
	.type	MatMult_SeqBAIJ_3, @function
MatMult_SeqBAIJ_3:
.LFB68:
	.loc 4 332 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$240, %rsp
	movq	%rdi, -200(%rbp)
	movq	%rsi, -208(%rbp)
	movq	%rdx, -216(%rbp)
	.loc 4 333 0
	movq	-200(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -168(%rbp)
	.loc 4 334 0
	movq	$0, -160(%rbp)
	.loc 4 338 0
	movq	$0, -48(%rbp)
	movl	$0, -40(%rbp)
	.loc 4 339 0
	movq	-168(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 346 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L302
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L303
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.25809, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$346, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L301
.L302:
	nop
	jmp	.L301
.L303:
	nop
.L301:
	.loc 4 347 0
	leaq	-184(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L281
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25809, %edx
	movl	$347, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L282
.L281:
	.loc 4 348 0
	leaq	-176(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L283
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25809, %edx
	movl	$348, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L282
.L283:
	.loc 4 350 0
	movq	-168(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 351 0
	movq	-168(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 352 0
	cmpl	$0, -36(%rbp)
	je	.L284
	.loc 4 353 0
	movq	-168(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 354 0
	movq	-168(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 355 0
	movq	-168(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L285
.L284:
	.loc 4 357 0
	movq	-168(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 358 0
	movq	-168(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 359 0
	movq	-176(%rbp), %rax
	movq	%rax, -160(%rbp)
.L285:
	.loc 4 362 0
	movl	$0, -76(%rbp)
	jmp	.L286
.L295:
	.loc 4 363 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 364 0
	movl	$0, %eax
	movq	%rax, -152(%rbp)
	movl	$0, %eax
	movq	%rax, -144(%rbp)
	movl	$0, %eax
	movq	%rax, -136(%rbp)
	.loc 4 365 0
	cmpl	$0, -52(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -40(%rbp)
.LBB7:
	.loc 4 366 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L287
.L288:
	addq	$64, -32(%rbp)
.L287:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L288
.LBE7:
.LBB8:
	.loc 4 367 0
	movq	-96(%rbp), %rcx
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L289
.L290:
	addq	$64, -16(%rbp)
.L289:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L290
.LBE8:
	.loc 4 368 0
	movl	$0, -56(%rbp)
	jmp	.L291
.L292:
	.loc 4 369 0
	movq	-184(%rbp), %rcx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -104(%rbp)
	addq	$4, -72(%rbp)
	movq	-104(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-104(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-104(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 370 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 371 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 372 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 373 0
	addq	$72, -96(%rbp)
	.loc 4 368 0
	addl	$1, -56(%rbp)
.L291:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L292
	.loc 4 375 0
	cmpl	$0, -36(%rbp)
	je	.L293
	movq	-176(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -160(%rbp)
.L293:
	.loc 4 376 0
	movq	-160(%rbp), %rax
	movq	-152(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-160(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-144(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-160(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-136(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 377 0
	cmpl	$0, -36(%rbp)
	jne	.L294
	addq	$24, -160(%rbp)
.L294:
	.loc 4 362 0
	addl	$1, -76(%rbp)
.L286:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L295
	.loc 4 379 0
	leaq	-184(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L296
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25809, %edx
	movl	$379, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L282
.L296:
	.loc 4 380 0
	leaq	-176(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L297
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25809, %edx
	movl	$380, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L282
.L297:
	.loc 4 381 0
	movq	-168(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC16(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-40(%rbp), %xmm0, %xmm0
	vmovsd	.LC17(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L298
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.25809, %edx
	movl	$381, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L282
.L298:
	.loc 4 382 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L299
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L299
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L299:
	movl	$0, %eax
.L282:
	.loc 4 383 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE68:
	.size	MatMult_SeqBAIJ_3, .-MatMult_SeqBAIJ_3
.globl MatMult_SeqBAIJ_4
	.type	MatMult_SeqBAIJ_4, @function
MatMult_SeqBAIJ_4:
.LFB69:
	.loc 4 388 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$256, %rsp
	movq	%rdi, -216(%rbp)
	movq	%rsi, -224(%rbp)
	movq	%rdx, -232(%rbp)
	.loc 4 389 0
	movq	-216(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -184(%rbp)
	.loc 4 390 0
	movq	$0, -176(%rbp)
	.loc 4 394 0
	movq	$0, -48(%rbp)
	movl	$0, -40(%rbp)
	.loc 4 395 0
	movq	-184(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 397 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L327
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L328
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.26011, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$397, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L326
.L327:
	nop
	jmp	.L326
.L328:
	nop
.L326:
	.loc 4 398 0
	leaq	-200(%rbp), %rdx
	movq	-224(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L306
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26011, %edx
	movl	$398, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L307
.L306:
	.loc 4 399 0
	leaq	-192(%rbp), %rdx
	movq	-232(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L308
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26011, %edx
	movl	$399, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L307
.L308:
	.loc 4 401 0
	movq	-184(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 402 0
	movq	-184(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 403 0
	cmpl	$0, -36(%rbp)
	je	.L309
	.loc 4 404 0
	movq	-184(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 405 0
	movq	-184(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 406 0
	movq	-184(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L310
.L309:
	.loc 4 408 0
	movq	-184(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 409 0
	movq	-184(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 410 0
	movq	-192(%rbp), %rax
	movq	%rax, -176(%rbp)
.L310:
	.loc 4 413 0
	movl	$0, -76(%rbp)
	jmp	.L311
.L320:
	.loc 4 414 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 415 0
	movl	$0, %eax
	movq	%rax, -168(%rbp)
	movl	$0, %eax
	movq	%rax, -160(%rbp)
	movl	$0, %eax
	movq	%rax, -152(%rbp)
	movl	$0, %eax
	movq	%rax, -144(%rbp)
	.loc 4 416 0
	cmpl	$0, -52(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -40(%rbp)
.LBB9:
	.loc 4 417 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L312
.L313:
	addq	$64, -32(%rbp)
.L312:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L313
.LBE9:
.LBB10:
	.loc 4 418 0
	movq	-96(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$7, %rdx
	addq	%rdx, %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$8, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L314
.L315:
	addq	$64, -16(%rbp)
.L314:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L315
.LBE10:
	.loc 4 419 0
	movl	$0, -56(%rbp)
	jmp	.L316
.L317:
	.loc 4 420 0
	movq	-200(%rbp), %rdx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$5, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -104(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 421 0
	movq	-104(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-104(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-104(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-104(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 422 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 423 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 424 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 425 0
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 426 0
	subq	$-128, -96(%rbp)
	.loc 4 419 0
	addl	$1, -56(%rbp)
.L316:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L317
	.loc 4 428 0
	cmpl	$0, -36(%rbp)
	je	.L318
	movq	-192(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$5, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -176(%rbp)
.L318:
	.loc 4 429 0
	movq	-176(%rbp), %rax
	movq	-168(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-176(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-160(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-176(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-152(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-176(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-144(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 430 0
	cmpl	$0, -36(%rbp)
	jne	.L319
	addq	$32, -176(%rbp)
.L319:
	.loc 4 413 0
	addl	$1, -76(%rbp)
.L311:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L320
	.loc 4 432 0
	leaq	-200(%rbp), %rdx
	movq	-224(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L321
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26011, %edx
	movl	$432, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L307
.L321:
	.loc 4 433 0
	leaq	-192(%rbp), %rdx
	movq	-232(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L322
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26011, %edx
	movl	$433, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L307
.L322:
	.loc 4 434 0
	movq	-184(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC18(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-40(%rbp), %xmm0, %xmm0
	vmovsd	.LC19(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L323
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26011, %edx
	movl	$434, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L307
.L323:
	.loc 4 435 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L324
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L324
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L324:
	movl	$0, %eax
.L307:
	.loc 4 436 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE69:
	.size	MatMult_SeqBAIJ_4, .-MatMult_SeqBAIJ_4
.globl MatMult_SeqBAIJ_5
	.type	MatMult_SeqBAIJ_5, @function
MatMult_SeqBAIJ_5:
.LFB70:
	.loc 4 441 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$272, %rsp
	movq	%rdi, -232(%rbp)
	movq	%rsi, -240(%rbp)
	movq	%rdx, -248(%rbp)
	.loc 4 442 0
	movq	-232(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -200(%rbp)
	.loc 4 443 0
	movq	$0, -112(%rbp)
	.loc 4 447 0
	movq	$0, -64(%rbp)
	.loc 4 448 0
	movl	$0, -40(%rbp)
	.loc 4 449 0
	movq	-200(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 451 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L352
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L353
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.26244, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$451, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L351
.L352:
	nop
	jmp	.L351
.L353:
	nop
.L351:
	.loc 4 452 0
	leaq	-216(%rbp), %rdx
	movq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L331
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26244, %edx
	movl	$452, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L332
.L331:
	.loc 4 453 0
	leaq	-208(%rbp), %rdx
	movq	-248(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L333
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26244, %edx
	movl	$453, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L332
.L333:
	.loc 4 455 0
	movq	-200(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 456 0
	movq	-200(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 457 0
	cmpl	$0, -36(%rbp)
	je	.L334
	.loc 4 458 0
	movq	-200(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -56(%rbp)
	.loc 4 459 0
	movq	-200(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 460 0
	movq	-200(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -64(%rbp)
	jmp	.L335
.L334:
	.loc 4 462 0
	movq	-200(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -56(%rbp)
	.loc 4 463 0
	movq	-200(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 464 0
	movq	-208(%rbp), %rax
	movq	%rax, -112(%rbp)
.L335:
	.loc 4 467 0
	movl	$0, -52(%rbp)
	jmp	.L336
.L345:
	.loc 4 468 0
	movq	-72(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -44(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 469 0
	movl	$0, %eax
	movq	%rax, -192(%rbp)
	movl	$0, %eax
	movq	%rax, -184(%rbp)
	movl	$0, %eax
	movq	%rax, -176(%rbp)
	movl	$0, %eax
	movq	%rax, -168(%rbp)
	movl	$0, %eax
	movq	%rax, -160(%rbp)
	.loc 4 470 0
	cmpl	$0, -44(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -40(%rbp)
.LBB11:
	.loc 4 471 0
	movq	-80(%rbp), %rax
	movl	-44(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-44(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-80(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L337
.L338:
	addq	$64, -32(%rbp)
.L337:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L338
.LBE11:
.LBB12:
	.loc 4 472 0
	movq	-96(%rbp), %rcx
	movl	-44(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	leaq	0(,%rax,4), %rdx
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -16(%rbp)
	movl	-44(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	leaq	0(,%rax,4), %rdx
	addq	%rdx, %rax
	salq	$4, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L339
.L340:
	addq	$64, -16(%rbp)
.L339:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L340
.LBE12:
	.loc 4 473 0
	movl	$0, -48(%rbp)
	jmp	.L341
.L342:
	.loc 4 474 0
	movq	-216(%rbp), %rcx
	movq	-80(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -104(%rbp)
	addq	$4, -80(%rbp)
	.loc 4 475 0
	movq	-104(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-104(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-104(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-104(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-104(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 476 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 477 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 478 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 479 0
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 480 0
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 481 0
	addq	$200, -96(%rbp)
	.loc 4 473 0
	addl	$1, -48(%rbp)
.L341:
	movl	-48(%rbp), %eax
	cmpl	-44(%rbp), %eax
	jl	.L342
	.loc 4 483 0
	cmpl	$0, -36(%rbp)
	je	.L343
	movq	-208(%rbp), %rcx
	movl	-52(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -112(%rbp)
.L343:
	.loc 4 484 0
	movq	-112(%rbp), %rax
	movq	-192(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-112(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-112(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-112(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-168(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-112(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-160(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 485 0
	cmpl	$0, -36(%rbp)
	jne	.L344
	addq	$40, -112(%rbp)
.L344:
	.loc 4 467 0
	addl	$1, -52(%rbp)
.L336:
	movl	-52(%rbp), %eax
	cmpl	-56(%rbp), %eax
	jl	.L345
	.loc 4 487 0
	leaq	-216(%rbp), %rdx
	movq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L346
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26244, %edx
	movl	$487, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L332
.L346:
	.loc 4 488 0
	leaq	-208(%rbp), %rdx
	movq	-248(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L347
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26244, %edx
	movl	$488, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L332
.L347:
	.loc 4 489 0
	movq	-200(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC20(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-40(%rbp), %xmm0, %xmm0
	vmovsd	.LC21(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L348
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26244, %edx
	movl	$489, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L332
.L348:
	.loc 4 490 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L349
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L349
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L349:
	movl	$0, %eax
.L332:
	.loc 4 491 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE70:
	.size	MatMult_SeqBAIJ_5, .-MatMult_SeqBAIJ_5
.globl MatMult_SeqBAIJ_6
	.type	MatMult_SeqBAIJ_6, @function
MatMult_SeqBAIJ_6:
.LFB71:
	.loc 4 497 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$288, %rsp
	movq	%rdi, -248(%rbp)
	movq	%rsi, -256(%rbp)
	movq	%rdx, -264(%rbp)
	.loc 4 498 0
	movq	-248(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -216(%rbp)
	.loc 4 499 0
	movq	$0, -208(%rbp)
	.loc 4 504 0
	movq	-216(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	movq	$0, -48(%rbp)
	movl	$0, -40(%rbp)
	.loc 4 505 0
	movq	-216(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 507 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L377
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L378
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.26516, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$507, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L376
.L377:
	nop
	jmp	.L376
.L378:
	nop
.L376:
	.loc 4 508 0
	leaq	-224(%rbp), %rdx
	movq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L356
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26516, %edx
	movl	$508, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L357
.L356:
	.loc 4 509 0
	leaq	-232(%rbp), %rdx
	movq	-264(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L358
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26516, %edx
	movl	$509, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L357
.L358:
	.loc 4 511 0
	movq	-216(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 512 0
	movq	-216(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 513 0
	cmpl	$0, -36(%rbp)
	je	.L359
	.loc 4 514 0
	movq	-216(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 515 0
	movq	-216(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 516 0
	movq	-216(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L360
.L359:
	.loc 4 518 0
	movq	-216(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 519 0
	movq	-216(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 520 0
	movq	-232(%rbp), %rax
	movq	%rax, -208(%rbp)
.L360:
	.loc 4 523 0
	movl	$0, -76(%rbp)
	jmp	.L361
.L370:
	.loc 4 524 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 525 0
	movl	$0, %eax
	movq	%rax, -200(%rbp)
	movl	$0, %eax
	movq	%rax, -192(%rbp)
	movl	$0, %eax
	movq	%rax, -184(%rbp)
	movl	$0, %eax
	movq	%rax, -176(%rbp)
	movl	$0, %eax
	movq	%rax, -168(%rbp)
	movl	$0, %eax
	movq	%rax, -160(%rbp)
	.loc 4 526 0
	cmpl	$0, -52(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -40(%rbp)
.LBB13:
	.loc 4 527 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L362
.L363:
	addq	$64, -32(%rbp)
.L362:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L363
.LBE13:
.LBB14:
	.loc 4 528 0
	movq	-96(%rbp), %rcx
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$5, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$6, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L364
.L365:
	addq	$64, -16(%rbp)
.L364:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L365
.LBE14:
	.loc 4 529 0
	movl	$0, -56(%rbp)
	jmp	.L366
.L367:
	.loc 4 530 0
	movq	-224(%rbp), %rcx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -152(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 531 0
	movq	-152(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-152(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-152(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-152(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-152(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-152(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 532 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 533 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 534 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 535 0
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 536 0
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 537 0
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 538 0
	addq	$288, -96(%rbp)
	.loc 4 529 0
	addl	$1, -56(%rbp)
.L366:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L367
	.loc 4 540 0
	cmpl	$0, -36(%rbp)
	je	.L368
	movq	-232(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -208(%rbp)
.L368:
	.loc 4 541 0
	movq	-208(%rbp), %rax
	movq	-200(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-208(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-192(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-168(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-160(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 542 0
	cmpl	$0, -36(%rbp)
	jne	.L369
	addq	$48, -208(%rbp)
.L369:
	.loc 4 523 0
	addl	$1, -76(%rbp)
.L361:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L370
	.loc 4 545 0
	leaq	-224(%rbp), %rdx
	movq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L371
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26516, %edx
	movl	$545, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L357
.L371:
	.loc 4 546 0
	leaq	-232(%rbp), %rdx
	movq	-264(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L372
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26516, %edx
	movl	$546, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L357
.L372:
	.loc 4 547 0
	movq	-216(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC22(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-40(%rbp), %xmm0, %xmm0
	vmovsd	.LC23(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L373
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26516, %edx
	movl	$547, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L357
.L373:
	.loc 4 548 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L374
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L374
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L374:
	movl	$0, %eax
.L357:
	.loc 4 549 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE71:
	.size	MatMult_SeqBAIJ_6, .-MatMult_SeqBAIJ_6
.globl MatMult_SeqBAIJ_7
	.type	MatMult_SeqBAIJ_7, @function
MatMult_SeqBAIJ_7:
.LFB72:
	.loc 4 554 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$312, %rsp
	movq	%rdi, -280(%rbp)
	movq	%rsi, -288(%rbp)
	movq	%rdx, -296(%rbp)
	.loc 4 555 0
	movq	-280(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -248(%rbp)
	.loc 4 556 0
	movq	$0, -240(%rbp)
	.loc 4 561 0
	movq	$0, -64(%rbp)
	movl	$0, -56(%rbp)
	.loc 4 562 0
	movq	-248(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -52(%rbp)
	.loc 4 564 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L402
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L403
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.26835, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$564, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L401
.L402:
	nop
	jmp	.L401
.L403:
	nop
.L401:
	.loc 4 565 0
	leaq	-256(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L381
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26835, %edx
	movl	$565, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L382
.L381:
	.loc 4 566 0
	leaq	-264(%rbp), %rdx
	movq	-296(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L383
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26835, %edx
	movl	$566, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L382
.L383:
	.loc 4 568 0
	movq	-248(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -88(%rbp)
	.loc 4 569 0
	movq	-248(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 570 0
	cmpl	$0, -52(%rbp)
	je	.L384
	.loc 4 571 0
	movq	-248(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -96(%rbp)
	.loc 4 572 0
	movq	-248(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 573 0
	movq	-248(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -64(%rbp)
	jmp	.L385
.L384:
	.loc 4 575 0
	movq	-248(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -96(%rbp)
	.loc 4 576 0
	movq	-248(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 577 0
	movq	-264(%rbp), %rax
	movq	%rax, -240(%rbp)
.L385:
	.loc 4 580 0
	movl	$0, -92(%rbp)
	jmp	.L386
.L395:
	.loc 4 581 0
	movq	-80(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-80(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -68(%rbp)
	addq	$4, -80(%rbp)
	.loc 4 582 0
	movl	$0, %eax
	movq	%rax, -232(%rbp)
	movl	$0, %eax
	movq	%rax, -224(%rbp)
	movl	$0, %eax
	movq	%rax, -216(%rbp)
	movl	$0, %eax
	movq	%rax, -208(%rbp)
	movl	$0, %eax
	movq	%rax, -200(%rbp)
	movl	$0, %eax
	movq	%rax, -192(%rbp)
	movl	$0, %eax
	movq	%rax, -184(%rbp)
	.loc 4 583 0
	cmpl	$0, -68(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -56(%rbp)
.LBB15:
	.loc 4 584 0
	movq	-88(%rbp), %rax
	movl	-68(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -48(%rbp)
	movl	-68(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-88(%rbp), %rax
	movq	%rax, -40(%rbp)
	jmp	.L387
.L388:
	addq	$64, -48(%rbp)
.L387:
	movq	-48(%rbp), %rax
	cmpq	-40(%rbp), %rax
	jb	.L388
.LBE15:
.LBB16:
	.loc 4 585 0
	movq	-112(%rbp), %rdx
	movl	-68(%rbp), %eax
	cltq
	imulq	$392, %rax, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -32(%rbp)
	movl	-68(%rbp), %eax
	cltq
	imulq	$784, %rax, %rax
	addq	-112(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L389
.L390:
	addq	$64, -32(%rbp)
.L389:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L390
.LBE16:
	.loc 4 586 0
	movl	$0, -72(%rbp)
	jmp	.L391
.L392:
	.loc 4 587 0
	movq	-256(%rbp), %rdx
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	0(,%rax,8), %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -176(%rbp)
	addq	$4, -88(%rbp)
	.loc 4 588 0
	movq	-176(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-176(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-176(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-176(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-176(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-176(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-176(%rbp), %rax
	addq	$48, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 589 0
	movq	-112(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 590 0
	movq	-112(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 591 0
	movq	-112(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 592 0
	movq	-112(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 593 0
	movq	-112(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 594 0
	movq	-112(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 595 0
	movq	-112(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 596 0
	addq	$392, -112(%rbp)
	.loc 4 586 0
	addl	$1, -72(%rbp)
.L391:
	movl	-72(%rbp), %eax
	cmpl	-68(%rbp), %eax
	jl	.L392
	.loc 4 598 0
	cmpl	$0, -52(%rbp)
	je	.L393
	movq	-264(%rbp), %rdx
	movl	-92(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	0(,%rax,8), %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -240(%rbp)
.L393:
	.loc 4 599 0
	movq	-240(%rbp), %rax
	movq	-232(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-240(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-224(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-216(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-208(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-192(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	48(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 600 0
	cmpl	$0, -52(%rbp)
	jne	.L394
	addq	$56, -240(%rbp)
.L394:
	.loc 4 580 0
	addl	$1, -92(%rbp)
.L386:
	movl	-92(%rbp), %eax
	cmpl	-96(%rbp), %eax
	jl	.L395
	.loc 4 603 0
	leaq	-256(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L396
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26835, %edx
	movl	$603, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L382
.L396:
	.loc 4 604 0
	leaq	-264(%rbp), %rdx
	movq	-296(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L397
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26835, %edx
	movl	$604, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L382
.L397:
	.loc 4 605 0
	movq	-248(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC24(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-56(%rbp), %xmm0, %xmm0
	vmovsd	.LC25(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L398
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.26835, %edx
	movl	$605, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L382
.L398:
	.loc 4 606 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L399
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L399
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L399:
	movl	$0, %eax
.L382:
	.loc 4 607 0
	addq	$312, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE72:
	.size	MatMult_SeqBAIJ_7, .-MatMult_SeqBAIJ_7
.globl MatMult_SeqBAIJ_15_ver1
	.type	MatMult_SeqBAIJ_15_ver1, @function
MatMult_SeqBAIJ_15_ver1:
.LFB73:
	.loc 4 615 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$312, %rsp
	movq	%rdi, -280(%rbp)
	movq	%rsi, -288(%rbp)
	movq	%rdx, -296(%rbp)
	.loc 4 616 0
	movq	-280(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -248(%rbp)
	.loc 4 617 0
	movq	$0, -240(%rbp)
	.loc 4 622 0
	movq	-248(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 623 0
	movq	$0, -32(%rbp)
	movl	$0, -24(%rbp)
	.loc 4 624 0
	movq	-248(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 626 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L425
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L426
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.27211, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$626, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L424
.L425:
	nop
	jmp	.L424
.L426:
	nop
.L424:
	.loc 4 627 0
	leaq	-256(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L406
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27211, %edx
	movl	$627, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L407
.L406:
	.loc 4 628 0
	leaq	-264(%rbp), %rdx
	movq	-296(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L408
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27211, %edx
	movl	$628, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L407
.L408:
	.loc 4 630 0
	movq	-248(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 631 0
	cmpl	$0, -20(%rbp)
	je	.L409
	.loc 4 632 0
	movq	-248(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -52(%rbp)
	.loc 4 633 0
	movq	-248(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 634 0
	movq	-248(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -32(%rbp)
	jmp	.L410
.L409:
	.loc 4 636 0
	movq	-248(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -52(%rbp)
	.loc 4 637 0
	movq	-248(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 638 0
	movq	-264(%rbp), %rax
	movq	%rax, -240(%rbp)
.L410:
	.loc 4 641 0
	movl	$0, -48(%rbp)
	jmp	.L411
.L418:
	.loc 4 642 0
	movl	-48(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-80(%rbp), %rax
	movl	(%rax), %edx
	movl	-48(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-80(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -36(%rbp)
	.loc 4 643 0
	movl	-48(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-80(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -64(%rbp)
	.loc 4 644 0
	movl	$0, %eax
	movq	%rax, -232(%rbp)
	movl	$0, %eax
	movq	%rax, -224(%rbp)
	movl	$0, %eax
	movq	%rax, -216(%rbp)
	movl	$0, %eax
	movq	%rax, -208(%rbp)
	movl	$0, %eax
	movq	%rax, -200(%rbp)
	movl	$0, %eax
	movq	%rax, -192(%rbp)
	movl	$0, %eax
	movq	%rax, -184(%rbp)
	.loc 4 645 0
	movl	$0, %eax
	movq	%rax, -176(%rbp)
	movl	$0, %eax
	movq	%rax, -168(%rbp)
	movl	$0, %eax
	movq	%rax, -160(%rbp)
	movl	$0, %eax
	movq	%rax, -152(%rbp)
	movl	$0, %eax
	movq	%rax, -144(%rbp)
	movl	$0, %eax
	movq	%rax, -136(%rbp)
	movl	$0, %eax
	movq	%rax, -128(%rbp)
	movl	$0, %eax
	movq	%rax, -120(%rbp)
	.loc 4 647 0
	cmpl	$0, -36(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -24(%rbp)
	.loc 4 648 0
	movl	$0, -44(%rbp)
	jmp	.L412
.L415:
	.loc 4 649 0
	movq	-256(%rbp), %rdx
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 651 0
	movl	$0, -40(%rbp)
	jmp	.L413
.L414:
	.loc 4 652 0
	movl	-40(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-112(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 653 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 654 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 655 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 656 0
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 657 0
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 658 0
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 659 0
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 660 0
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 661 0
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 662 0
	movq	-96(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 663 0
	movq	-96(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 664 0
	movq	-96(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 665 0
	movq	-96(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 666 0
	movq	-96(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-128(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -128(%rbp)
	.loc 4 667 0
	movq	-96(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vmovsd	-120(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -120(%rbp)
	.loc 4 668 0
	addq	$120, -96(%rbp)
	.loc 4 651 0
	addl	$1, -40(%rbp)
.L413:
	cmpl	$14, -40(%rbp)
	jle	.L414
	.loc 4 648 0
	addl	$1, -44(%rbp)
.L412:
	movl	-44(%rbp), %eax
	cmpl	-36(%rbp), %eax
	jl	.L415
	.loc 4 671 0
	cmpl	$0, -20(%rbp)
	je	.L416
	movq	-264(%rbp), %rdx
	movl	-48(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -240(%rbp)
.L416:
	.loc 4 672 0
	movq	-240(%rbp), %rax
	movq	-232(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-240(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-224(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-216(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-208(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-192(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	48(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 673 0
	movq	-240(%rbp), %rax
	leaq	56(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	64(%rax), %rdx
	movq	-168(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	72(%rax), %rdx
	movq	-160(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	80(%rax), %rdx
	movq	-152(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	88(%rax), %rdx
	movq	-144(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	96(%rax), %rdx
	movq	-136(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	104(%rax), %rdx
	movq	-128(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	112(%rax), %rdx
	movq	-120(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 675 0
	cmpl	$0, -20(%rbp)
	jne	.L417
	addq	$120, -240(%rbp)
.L417:
	.loc 4 641 0
	addl	$1, -48(%rbp)
.L411:
	movl	-48(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L418
	.loc 4 678 0
	leaq	-256(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L419
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27211, %edx
	movl	$678, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L407
.L419:
	.loc 4 679 0
	leaq	-264(%rbp), %rdx
	movq	-296(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L420
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27211, %edx
	movl	$679, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L407
.L420:
	.loc 4 680 0
	movq	-248(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC26(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-24(%rbp), %xmm0, %xmm0
	vmovsd	.LC27(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L421
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27211, %edx
	movl	$680, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L407
.L421:
	.loc 4 681 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L422
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L422
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L422:
	movl	$0, %eax
.L407:
	.loc 4 682 0
	addq	$312, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE73:
	.size	MatMult_SeqBAIJ_15_ver1, .-MatMult_SeqBAIJ_15_ver1
.globl MatMult_SeqBAIJ_15_ver2
	.type	MatMult_SeqBAIJ_15_ver2, @function
MatMult_SeqBAIJ_15_ver2:
.LFB74:
	.loc 4 688 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$328, %rsp
	movq	%rdi, -296(%rbp)
	movq	%rsi, -304(%rbp)
	movq	%rdx, -312(%rbp)
	.loc 4 689 0
	movq	-296(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -264(%rbp)
	.loc 4 690 0
	movq	$0, -256(%rbp)
	.loc 4 695 0
	movq	-264(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 696 0
	movq	$0, -32(%rbp)
	movl	$0, -24(%rbp)
	.loc 4 697 0
	movq	-264(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 699 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L446
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L447
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.27446, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$699, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L445
.L446:
	nop
	jmp	.L445
.L447:
	nop
.L445:
	.loc 4 700 0
	leaq	-272(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L429
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27446, %edx
	movl	$700, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L430
.L429:
	.loc 4 701 0
	leaq	-280(%rbp), %rdx
	movq	-312(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L431
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27446, %edx
	movl	$701, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L430
.L431:
	.loc 4 703 0
	movq	-264(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -88(%rbp)
	.loc 4 704 0
	cmpl	$0, -20(%rbp)
	je	.L432
	.loc 4 705 0
	movq	-264(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -48(%rbp)
	.loc 4 706 0
	movq	-264(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 707 0
	movq	-264(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -32(%rbp)
	jmp	.L433
.L432:
	.loc 4 709 0
	movq	-264(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -48(%rbp)
	.loc 4 710 0
	movq	-264(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 711 0
	movq	-280(%rbp), %rax
	movq	%rax, -256(%rbp)
.L433:
	.loc 4 714 0
	movl	$0, -44(%rbp)
	jmp	.L434
.L439:
	.loc 4 715 0
	movl	-44(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %edx
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -36(%rbp)
	.loc 4 716 0
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 717 0
	movl	$0, %eax
	movq	%rax, -248(%rbp)
	movl	$0, %eax
	movq	%rax, -240(%rbp)
	movl	$0, %eax
	movq	%rax, -232(%rbp)
	movl	$0, %eax
	movq	%rax, -224(%rbp)
	movl	$0, %eax
	movq	%rax, -216(%rbp)
	movl	$0, %eax
	movq	%rax, -208(%rbp)
	movl	$0, %eax
	movq	%rax, -200(%rbp)
	.loc 4 718 0
	movl	$0, %eax
	movq	%rax, -192(%rbp)
	movl	$0, %eax
	movq	%rax, -184(%rbp)
	movl	$0, %eax
	movq	%rax, -176(%rbp)
	movl	$0, %eax
	movq	%rax, -168(%rbp)
	movl	$0, %eax
	movq	%rax, -160(%rbp)
	movl	$0, %eax
	movq	%rax, -152(%rbp)
	movl	$0, %eax
	movq	%rax, -144(%rbp)
	movl	$0, %eax
	movq	%rax, -136(%rbp)
	.loc 4 720 0
	cmpl	$0, -36(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -24(%rbp)
	.loc 4 721 0
	movl	$0, -40(%rbp)
	jmp	.L435
.L436:
	.loc 4 722 0
	movq	-272(%rbp), %rdx
	movl	-40(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -128(%rbp)
	.loc 4 723 0
	movq	-128(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-128(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-128(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	movq	-128(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 725 0
	movq	-88(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-248(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)
	.loc 4 726 0
	movq	-88(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-240(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -240(%rbp)
	.loc 4 727 0
	movq	-88(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 728 0
	movq	-88(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 729 0
	movq	-88(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$392, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 730 0
	movq	-88(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$400, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 731 0
	movq	-88(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$408, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 732 0
	movq	-88(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$416, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 733 0
	movq	-88(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$424, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 734 0
	movq	-88(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$432, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 735 0
	movq	-88(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$440, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 736 0
	movq	-88(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$448, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 737 0
	movq	-88(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$456, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 738 0
	movq	-88(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$464, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 739 0
	movq	-88(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$472, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 741 0
	addq	$480, -88(%rbp)
	.loc 4 743 0
	movq	-128(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-128(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-128(%rbp), %rax
	addq	$48, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	movq	-128(%rbp), %rax
	addq	$56, %rax
	movq	(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 745 0
	movq	-88(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-248(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)
	.loc 4 746 0
	movq	-88(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-240(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -240(%rbp)
	.loc 4 747 0
	movq	-88(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 748 0
	movq	-88(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 749 0
	movq	-88(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$392, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 750 0
	movq	-88(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$400, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 751 0
	movq	-88(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$408, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 752 0
	movq	-88(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$416, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 753 0
	movq	-88(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$424, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 754 0
	movq	-88(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$432, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 755 0
	movq	-88(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$440, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 756 0
	movq	-88(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$448, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 757 0
	movq	-88(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$456, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 758 0
	movq	-88(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$464, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 759 0
	movq	-88(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$472, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 760 0
	addq	$480, -88(%rbp)
	.loc 4 762 0
	movq	-128(%rbp), %rax
	addq	$64, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-128(%rbp), %rax
	addq	$72, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-128(%rbp), %rax
	addq	$80, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	movq	-128(%rbp), %rax
	addq	$88, %rax
	movq	(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 763 0
	movq	-88(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-248(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)
	.loc 4 764 0
	movq	-88(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-240(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -240(%rbp)
	.loc 4 765 0
	movq	-88(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 766 0
	movq	-88(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 767 0
	movq	-88(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$392, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 768 0
	movq	-88(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$400, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 769 0
	movq	-88(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$408, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 770 0
	movq	-88(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$416, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 771 0
	movq	-88(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$424, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 772 0
	movq	-88(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$432, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 773 0
	movq	-88(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$440, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 774 0
	movq	-88(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$448, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 775 0
	movq	-88(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$456, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 776 0
	movq	-88(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$464, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 777 0
	movq	-88(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$472, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 778 0
	addq	$480, -88(%rbp)
	.loc 4 780 0
	movq	-128(%rbp), %rax
	addq	$96, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-128(%rbp), %rax
	addq	$104, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-128(%rbp), %rax
	addq	$112, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 781 0
	movq	-88(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-248(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)
	.loc 4 782 0
	movq	-88(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-240(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -240(%rbp)
	.loc 4 783 0
	movq	-88(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 784 0
	movq	-88(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 785 0
	movq	-88(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 786 0
	movq	-88(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 787 0
	movq	-88(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 788 0
	movq	-88(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 789 0
	movq	-88(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 790 0
	movq	-88(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 791 0
	movq	-88(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 792 0
	movq	-88(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 793 0
	movq	-88(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 794 0
	movq	-88(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 795 0
	movq	-88(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 796 0
	addq	$360, -88(%rbp)
	.loc 4 721 0
	addl	$1, -40(%rbp)
.L435:
	movl	-40(%rbp), %eax
	cmpl	-36(%rbp), %eax
	jl	.L436
	.loc 4 798 0
	cmpl	$0, -20(%rbp)
	je	.L437
	movq	-280(%rbp), %rdx
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -256(%rbp)
.L437:
	.loc 4 799 0
	movq	-256(%rbp), %rax
	movq	-248(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-256(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-240(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-232(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-224(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-216(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-208(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	48(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 800 0
	movq	-256(%rbp), %rax
	leaq	56(%rax), %rdx
	movq	-192(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	64(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	72(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	80(%rax), %rdx
	movq	-168(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	88(%rax), %rdx
	movq	-160(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	96(%rax), %rdx
	movq	-152(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	104(%rax), %rdx
	movq	-144(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-256(%rbp), %rax
	leaq	112(%rax), %rdx
	movq	-136(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 802 0
	cmpl	$0, -20(%rbp)
	jne	.L438
	addq	$120, -256(%rbp)
.L438:
	.loc 4 714 0
	addl	$1, -44(%rbp)
.L434:
	movl	-44(%rbp), %eax
	cmpl	-48(%rbp), %eax
	jl	.L439
	.loc 4 805 0
	leaq	-272(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L440
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27446, %edx
	movl	$805, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L430
.L440:
	.loc 4 806 0
	leaq	-280(%rbp), %rdx
	movq	-312(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L441
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27446, %edx
	movl	$806, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L430
.L441:
	.loc 4 807 0
	movq	-264(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC26(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-24(%rbp), %xmm0, %xmm0
	vmovsd	.LC27(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L442
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.27446, %edx
	movl	$807, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L430
.L442:
	.loc 4 808 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L443
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L443
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L443:
	movl	$0, %eax
.L430:
	.loc 4 809 0
	addq	$328, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE74:
	.size	MatMult_SeqBAIJ_15_ver2, .-MatMult_SeqBAIJ_15_ver2
.globl MatMult_SeqBAIJ_15_ver3
	.type	MatMult_SeqBAIJ_15_ver3, @function
MatMult_SeqBAIJ_15_ver3:
.LFB75:
	.loc 4 815 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$360, %rsp
	movq	%rdi, -328(%rbp)
	movq	%rsi, -336(%rbp)
	movq	%rdx, -344(%rbp)
	.loc 4 816 0
	movq	-328(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -296(%rbp)
	.loc 4 817 0
	movq	$0, -288(%rbp)
	.loc 4 822 0
	movq	-296(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 823 0
	movq	$0, -32(%rbp)
	movl	$0, -24(%rbp)
	.loc 4 824 0
	movq	-296(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 826 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L467
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L468
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.28485, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$826, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L466
.L467:
	nop
	jmp	.L466
.L468:
	nop
.L466:
	.loc 4 827 0
	leaq	-304(%rbp), %rdx
	movq	-336(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L450
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.28485, %edx
	movl	$827, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L451
.L450:
	.loc 4 828 0
	leaq	-312(%rbp), %rdx
	movq	-344(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L452
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.28485, %edx
	movl	$828, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L451
.L452:
	.loc 4 830 0
	movq	-296(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -88(%rbp)
	.loc 4 831 0
	cmpl	$0, -20(%rbp)
	je	.L453
	.loc 4 832 0
	movq	-296(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -48(%rbp)
	.loc 4 833 0
	movq	-296(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 834 0
	movq	-296(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -32(%rbp)
	jmp	.L454
.L453:
	.loc 4 836 0
	movq	-296(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -48(%rbp)
	.loc 4 837 0
	movq	-296(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 838 0
	movq	-312(%rbp), %rax
	movq	%rax, -288(%rbp)
.L454:
	.loc 4 841 0
	movl	$0, -44(%rbp)
	jmp	.L455
.L460:
	.loc 4 842 0
	movl	-44(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %edx
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -36(%rbp)
	.loc 4 843 0
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 844 0
	movl	$0, %eax
	movq	%rax, -280(%rbp)
	movl	$0, %eax
	movq	%rax, -272(%rbp)
	movl	$0, %eax
	movq	%rax, -264(%rbp)
	movl	$0, %eax
	movq	%rax, -256(%rbp)
	movl	$0, %eax
	movq	%rax, -248(%rbp)
	movl	$0, %eax
	movq	%rax, -240(%rbp)
	movl	$0, %eax
	movq	%rax, -232(%rbp)
	.loc 4 845 0
	movl	$0, %eax
	movq	%rax, -224(%rbp)
	movl	$0, %eax
	movq	%rax, -216(%rbp)
	movl	$0, %eax
	movq	%rax, -208(%rbp)
	movl	$0, %eax
	movq	%rax, -200(%rbp)
	movl	$0, %eax
	movq	%rax, -192(%rbp)
	movl	$0, %eax
	movq	%rax, -184(%rbp)
	movl	$0, %eax
	movq	%rax, -176(%rbp)
	movl	$0, %eax
	movq	%rax, -168(%rbp)
	.loc 4 847 0
	cmpl	$0, -36(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -24(%rbp)
	.loc 4 848 0
	movl	$0, -40(%rbp)
	jmp	.L456
.L457:
	.loc 4 849 0
	movq	-304(%rbp), %rdx
	movl	-40(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -160(%rbp)
	.loc 4 850 0
	movq	-160(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-160(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-160(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-160(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-160(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-160(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-160(%rbp), %rax
	addq	$48, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 851 0
	movq	-160(%rbp), %rax
	addq	$56, %rax
	movq	(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 853 0
	movq	-88(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$480, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$600, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$720, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$840, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-280(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -280(%rbp)
	.loc 4 854 0
	movq	-88(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$488, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$608, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$728, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$848, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-272(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -272(%rbp)
	.loc 4 855 0
	movq	-88(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$496, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$616, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$736, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$856, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-264(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -264(%rbp)
	.loc 4 856 0
	movq	-88(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$504, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$624, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$744, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$864, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-256(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -256(%rbp)
	.loc 4 857 0
	movq	-88(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$392, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$512, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$632, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$752, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$872, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-248(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)
	.loc 4 858 0
	movq	-88(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$400, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$520, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$640, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$760, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$880, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-240(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -240(%rbp)
	.loc 4 859 0
	movq	-88(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$408, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$528, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$648, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$768, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$888, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 860 0
	movq	-88(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$416, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$536, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$656, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$776, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$896, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 861 0
	movq	-88(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$424, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$544, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$664, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$784, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$904, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 862 0
	movq	-88(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$432, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$552, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$672, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$792, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$912, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 863 0
	movq	-88(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$440, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$560, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$680, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$800, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$920, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 864 0
	movq	-88(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$448, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$568, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$688, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$808, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$928, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 865 0
	movq	-88(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$456, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$576, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$696, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$816, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$936, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 866 0
	movq	-88(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$464, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$584, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$704, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$824, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$944, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 867 0
	movq	-88(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$472, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$592, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$712, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$832, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$952, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 868 0
	addq	$960, -88(%rbp)
	.loc 4 870 0
	movq	-160(%rbp), %rax
	addq	$64, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-160(%rbp), %rax
	addq	$72, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-160(%rbp), %rax
	addq	$80, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-160(%rbp), %rax
	addq	$88, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-160(%rbp), %rax
	addq	$96, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-160(%rbp), %rax
	addq	$104, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-160(%rbp), %rax
	addq	$112, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 872 0
	movq	-88(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$480, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$600, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$720, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-280(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -280(%rbp)
	.loc 4 873 0
	movq	-88(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$488, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$608, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$728, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-272(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -272(%rbp)
	.loc 4 874 0
	movq	-88(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$496, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$616, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$736, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-264(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -264(%rbp)
	.loc 4 875 0
	movq	-88(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$504, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$624, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$744, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-256(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -256(%rbp)
	.loc 4 876 0
	movq	-88(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$392, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$512, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$632, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$752, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-248(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)
	.loc 4 877 0
	movq	-88(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$400, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$520, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$640, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$760, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-240(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -240(%rbp)
	.loc 4 878 0
	movq	-88(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$408, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$528, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$648, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$768, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 879 0
	movq	-88(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$416, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$536, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$656, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$776, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 880 0
	movq	-88(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$424, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$544, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$664, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$784, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 881 0
	movq	-88(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$432, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$552, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$672, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$792, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 882 0
	movq	-88(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$440, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$560, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$680, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$800, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 883 0
	movq	-88(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$448, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$568, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$688, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$808, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 884 0
	movq	-88(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$456, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$576, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$696, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$816, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 885 0
	movq	-88(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$464, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$584, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$704, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$824, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 886 0
	movq	-88(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$472, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$592, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$712, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$832, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 887 0
	addq	$840, -88(%rbp)
	.loc 4 848 0
	addl	$1, -40(%rbp)
.L456:
	movl	-40(%rbp), %eax
	cmpl	-36(%rbp), %eax
	jl	.L457
	.loc 4 889 0
	cmpl	$0, -20(%rbp)
	je	.L458
	movq	-312(%rbp), %rdx
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -288(%rbp)
.L458:
	.loc 4 890 0
	movq	-288(%rbp), %rax
	movq	-280(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-288(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-272(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-264(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-256(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-248(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-240(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	48(%rax), %rdx
	movq	-232(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 891 0
	movq	-288(%rbp), %rax
	leaq	56(%rax), %rdx
	movq	-224(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	64(%rax), %rdx
	movq	-216(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	72(%rax), %rdx
	movq	-208(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	80(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	88(%rax), %rdx
	movq	-192(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	96(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	104(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-288(%rbp), %rax
	leaq	112(%rax), %rdx
	movq	-168(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 893 0
	cmpl	$0, -20(%rbp)
	jne	.L459
	addq	$120, -288(%rbp)
.L459:
	.loc 4 841 0
	addl	$1, -44(%rbp)
.L455:
	movl	-44(%rbp), %eax
	cmpl	-48(%rbp), %eax
	jl	.L460
	.loc 4 896 0
	leaq	-304(%rbp), %rdx
	movq	-336(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L461
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.28485, %edx
	movl	$896, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L451
.L461:
	.loc 4 897 0
	leaq	-312(%rbp), %rdx
	movq	-344(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L462
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.28485, %edx
	movl	$897, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L451
.L462:
	.loc 4 898 0
	movq	-296(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC26(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-24(%rbp), %xmm0, %xmm0
	vmovsd	.LC27(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L463
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.28485, %edx
	movl	$898, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L451
.L463:
	.loc 4 899 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L464
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L464
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L464:
	movl	$0, %eax
.L451:
	.loc 4 900 0
	addq	$360, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE75:
	.size	MatMult_SeqBAIJ_15_ver3, .-MatMult_SeqBAIJ_15_ver3
.globl MatMult_SeqBAIJ_15_ver4
	.type	MatMult_SeqBAIJ_15_ver4, @function
MatMult_SeqBAIJ_15_ver4:
.LFB76:
	.loc 4 907 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$408, %rsp
	movq	%rdi, -376(%rbp)
	movq	%rsi, -384(%rbp)
	movq	%rdx, -392(%rbp)
	.loc 4 908 0
	movq	-376(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -352(%rbp)
	.loc 4 909 0
	movq	$0, -344(%rbp)
	.loc 4 914 0
	movq	-352(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 915 0
	movq	$0, -32(%rbp)
	movl	$0, -24(%rbp)
	.loc 4 916 0
	movq	-352(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 918 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L488
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L489
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.29563, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$918, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L487
.L488:
	nop
	jmp	.L487
.L489:
	nop
.L487:
	.loc 4 919 0
	leaq	-360(%rbp), %rdx
	movq	-384(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L471
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.29563, %edx
	movl	$919, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L472
.L471:
	.loc 4 920 0
	leaq	-368(%rbp), %rdx
	movq	-392(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L473
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.29563, %edx
	movl	$920, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L472
.L473:
	.loc 4 922 0
	movq	-352(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -88(%rbp)
	.loc 4 923 0
	cmpl	$0, -20(%rbp)
	je	.L474
	.loc 4 924 0
	movq	-352(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -48(%rbp)
	.loc 4 925 0
	movq	-352(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 926 0
	movq	-352(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -32(%rbp)
	jmp	.L475
.L474:
	.loc 4 928 0
	movq	-352(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -48(%rbp)
	.loc 4 929 0
	movq	-352(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 930 0
	movq	-368(%rbp), %rax
	movq	%rax, -344(%rbp)
.L475:
	.loc 4 933 0
	movl	$0, -44(%rbp)
	jmp	.L476
.L481:
	.loc 4 934 0
	movl	-44(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %edx
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -36(%rbp)
	.loc 4 935 0
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 936 0
	movl	$0, %eax
	movq	%rax, -336(%rbp)
	movl	$0, %eax
	movq	%rax, -328(%rbp)
	movl	$0, %eax
	movq	%rax, -320(%rbp)
	movl	$0, %eax
	movq	%rax, -312(%rbp)
	movl	$0, %eax
	movq	%rax, -304(%rbp)
	movl	$0, %eax
	movq	%rax, -296(%rbp)
	movl	$0, %eax
	movq	%rax, -288(%rbp)
	.loc 4 937 0
	movl	$0, %eax
	movq	%rax, -280(%rbp)
	movl	$0, %eax
	movq	%rax, -272(%rbp)
	movl	$0, %eax
	movq	%rax, -264(%rbp)
	movl	$0, %eax
	movq	%rax, -256(%rbp)
	movl	$0, %eax
	movq	%rax, -248(%rbp)
	movl	$0, %eax
	movq	%rax, -240(%rbp)
	movl	$0, %eax
	movq	%rax, -232(%rbp)
	movl	$0, %eax
	movq	%rax, -224(%rbp)
	.loc 4 939 0
	cmpl	$0, -36(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -24(%rbp)
	.loc 4 940 0
	movl	$0, -40(%rbp)
	jmp	.L477
.L478:
	.loc 4 941 0
	movq	-360(%rbp), %rdx
	movl	-40(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -216(%rbp)
	.loc 4 942 0
	movq	-216(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -208(%rbp)
	movq	-216(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -200(%rbp)
	movq	-216(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -192(%rbp)
	movq	-216(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -184(%rbp)
	movq	-216(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	movq	-216(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-216(%rbp), %rax
	addq	$48, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	.loc 4 943 0
	movq	-216(%rbp), %rax
	addq	$56, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-216(%rbp), %rax
	addq	$64, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-216(%rbp), %rax
	addq	$72, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-216(%rbp), %rax
	addq	$80, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-216(%rbp), %rax
	addq	$88, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-216(%rbp), %rax
	addq	$96, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-216(%rbp), %rax
	addq	$104, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	movq	-216(%rbp), %rax
	addq	$112, %rax
	movq	(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 945 0
	movq	-88(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$480, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$600, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$720, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$840, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$960, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1080, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1440, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1560, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1680, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-336(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -336(%rbp)
	.loc 4 946 0
	movq	-88(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$488, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$608, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$728, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$848, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$968, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1088, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1448, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1568, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1688, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-328(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -328(%rbp)
	.loc 4 947 0
	movq	-88(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$496, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$616, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$736, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$856, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$976, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1096, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1456, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1576, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1696, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-320(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -320(%rbp)
	.loc 4 948 0
	movq	-88(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$504, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$624, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$744, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$864, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$984, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1464, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1584, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1704, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-312(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -312(%rbp)
	.loc 4 949 0
	movq	-88(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$392, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$512, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$632, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$752, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$872, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$992, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1472, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1592, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1712, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-304(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -304(%rbp)
	.loc 4 950 0
	movq	-88(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$400, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$520, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$640, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$760, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$880, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1000, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1480, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1600, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1720, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-296(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -296(%rbp)
	.loc 4 951 0
	movq	-88(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$408, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$528, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$648, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$768, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$888, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1008, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1488, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1608, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1728, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-288(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -288(%rbp)
	.loc 4 952 0
	movq	-88(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$416, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$536, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$656, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$776, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$896, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1016, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1496, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1616, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1736, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-280(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -280(%rbp)
	.loc 4 953 0
	movq	-88(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$424, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$544, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$664, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$784, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$904, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1024, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1504, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1624, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1744, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-272(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -272(%rbp)
	.loc 4 954 0
	movq	-88(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$432, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$552, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$672, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$792, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$912, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1032, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1392, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1512, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1632, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1752, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-264(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -264(%rbp)
	.loc 4 955 0
	movq	-88(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$440, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$560, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$680, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$800, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$920, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1040, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1400, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1520, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1640, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1760, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-256(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -256(%rbp)
	.loc 4 956 0
	movq	-88(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$448, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$568, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$688, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$808, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$928, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1048, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1408, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1528, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1648, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1768, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-248(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)
	.loc 4 957 0
	movq	-88(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$456, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$576, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$696, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$816, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$936, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1056, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1416, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1536, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1656, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1776, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-240(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -240(%rbp)
	.loc 4 958 0
	movq	-88(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$464, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$584, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$704, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$824, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$944, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1064, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1424, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1544, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1664, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1784, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-232(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -232(%rbp)
	.loc 4 959 0
	movq	-88(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-208(%rbp), %xmm0, %xmm1
	movq	-88(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-200(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-192(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$472, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-184(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$592, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$712, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$832, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$952, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1072, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1432, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1552, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1672, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-88(%rbp), %rax
	addq	$1792, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-96(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 960 0
	addq	$1800, -88(%rbp)
	.loc 4 940 0
	addl	$1, -40(%rbp)
.L477:
	movl	-40(%rbp), %eax
	cmpl	-36(%rbp), %eax
	jl	.L478
	.loc 4 962 0
	cmpl	$0, -20(%rbp)
	je	.L479
	movq	-368(%rbp), %rdx
	movl	-44(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	movq	%rax, %rcx
	salq	$4, %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -344(%rbp)
.L479:
	.loc 4 963 0
	movq	-344(%rbp), %rax
	movq	-336(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-344(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-328(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-320(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-312(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-304(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-296(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	48(%rax), %rdx
	movq	-288(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 964 0
	movq	-344(%rbp), %rax
	leaq	56(%rax), %rdx
	movq	-280(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	64(%rax), %rdx
	movq	-272(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	72(%rax), %rdx
	movq	-264(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	80(%rax), %rdx
	movq	-256(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	88(%rax), %rdx
	movq	-248(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	96(%rax), %rdx
	movq	-240(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	104(%rax), %rdx
	movq	-232(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-344(%rbp), %rax
	leaq	112(%rax), %rdx
	movq	-224(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 966 0
	cmpl	$0, -20(%rbp)
	jne	.L480
	addq	$120, -344(%rbp)
.L480:
	.loc 4 933 0
	addl	$1, -44(%rbp)
.L476:
	movl	-44(%rbp), %eax
	cmpl	-48(%rbp), %eax
	jl	.L481
	.loc 4 969 0
	leaq	-360(%rbp), %rdx
	movq	-384(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L482
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.29563, %edx
	movl	$969, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L472
.L482:
	.loc 4 970 0
	leaq	-368(%rbp), %rdx
	movq	-392(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L483
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.29563, %edx
	movl	$970, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L472
.L483:
	.loc 4 971 0
	movq	-352(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC26(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm1
	vcvtsi2sd	-24(%rbp), %xmm0, %xmm0
	vmovsd	.LC27(%rip), %xmm2
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -76(%rbp)
	cmpl	$0, -76(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L484
	movl	-76(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.29563, %edx
	movl	$971, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L472
.L484:
	.loc 4 972 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L485
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L485
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L485:
	movl	$0, %eax
.L472:
	.loc 4 973 0
	addq	$408, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE76:
	.size	MatMult_SeqBAIJ_15_ver4, .-MatMult_SeqBAIJ_15_ver4
	.section	.rodata
.LC29:
	.string	"N"
	.text
.globl MatMult_SeqBAIJ_N
	.type	MatMult_SeqBAIJ_N, @function
MatMult_SeqBAIJ_N:
.LFB77:
	.loc 4 982 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$264, %rsp
	movq	%rdi, -200(%rbp)
	movq	%rsi, -208(%rbp)
	movq	%rdx, -216(%rbp)
	.loc 4 983 0
	movq	-200(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -136(%rbp)
	.loc 4 984 0
	movq	$0, -128(%rbp)
	.loc 4 987 0
	movq	-136(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	movq	-200(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -56(%rbp)
	movq	-136(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -44(%rbp)
	.loc 4 988 0
	movq	$0, -32(%rbp)
	movl	$0, -24(%rbp)
	.loc 4 989 0
	movq	-136(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 991 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L514
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L515
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.30632, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$991, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L513
.L514:
	nop
	jmp	.L513
.L515:
	nop
.L513:
	.loc 4 992 0
	leaq	-144(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L492
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30632, %edx
	movl	$992, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L493
.L492:
	.loc 4 993 0
	leaq	-152(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L494
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30632, %edx
	movl	$993, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L493
.L494:
	.loc 4 995 0
	movq	-136(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 996 0
	movq	-136(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 997 0
	cmpl	$0, -20(%rbp)
	je	.L495
	.loc 4 998 0
	movq	-136(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 999 0
	movq	-136(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1000 0
	movq	-136(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -32(%rbp)
	jmp	.L496
.L495:
	.loc 4 1002 0
	movq	-136(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1003 0
	movq	-136(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1004 0
	movq	-152(%rbp), %rax
	movq	%rax, -128(%rbp)
.L496:
	.loc 4 1007 0
	movq	-136(%rbp), %rax
	movq	240(%rax), %rax
	testq	%rax, %rax
	jne	.L497
	.loc 4 1008 0
	movq	-200(%rbp), %rax
	movq	456(%rax), %rax
	movl	4(%rax), %edx
	movq	-200(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	cmpl	%eax, %edx
	cmovge	%edx, %eax
	movl	%eax, -36(%rbp)
	.loc 4 1009 0
	movl	-36(%rbp), %eax
	addl	$1, %eax
	cltq
	salq	$3, %rax
	testq	%rax, %rax
	je	.L498
	movq	PetscTrMalloc(%rip), %rax
	movq	-136(%rbp), %rdx
	addq	$240, %rdx
	movl	-36(%rbp), %ecx
	addl	$1, %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,8), %rbx
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30632, %edx
	movl	$1009, %esi
	movq	%rbx, %rdi
	call	*%rax
	jmp	.L499
.L498:
	movq	-136(%rbp), %rax
	movq	$0, 240(%rax)
	movl	$0, %eax
.L499:
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L497
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30632, %edx
	movl	$1009, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L493
.L497:
	.loc 4 1011 0
	movq	-136(%rbp), %rax
	movq	240(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 1012 0
	movl	$0, -76(%rbp)
	jmp	.L500
.L507:
	.loc 4 1013 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -48(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 1014 0
	movl	-48(%rbp), %eax
	imull	-56(%rbp), %eax
	movl	%eax, -40(%rbp)
	.loc 4 1015 0
	movq	-112(%rbp), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1016 0
	cmpl	$0, -48(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -24(%rbp)
	.loc 4 1017 0
	movl	$0, -52(%rbp)
	jmp	.L501
.L504:
	.loc 4 1018 0
	movq	-144(%rbp), %rdx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	imull	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -120(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 1019 0
	movl	$0, -36(%rbp)
	jmp	.L502
.L503:
	movl	-36(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-104(%rbp), %rax
	movl	-36(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	-120(%rbp), %rdx
	movq	(%rdx), %rdx
	movq	%rdx, (%rax)
	addl	$1, -36(%rbp)
.L502:
	movl	-36(%rbp), %eax
	cmpl	-56(%rbp), %eax
	jl	.L503
	.loc 4 1020 0
	movl	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -104(%rbp)
	.loc 4 1017 0
	addl	$1, -52(%rbp)
.L501:
	movl	-52(%rbp), %eax
	cmpl	-48(%rbp), %eax
	jl	.L504
	.loc 4 1022 0
	cmpl	$0, -20(%rbp)
	je	.L505
	movq	-152(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	imull	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -128(%rbp)
.L505:
.LBB17:
	.loc 4 1023 0
	movabsq	$4607182418800017408, %rax
	movq	%rax, -160(%rbp)
	movl	$0, %eax
	movq	%rax, -168(%rbp)
	movl	$1, -172(%rbp)
	movl	-56(%rbp), %eax
	movl	%eax, -176(%rbp)
	movl	-40(%rbp), %eax
	movl	%eax, -180(%rbp)
	leaq	-176(%rbp), %rdi
	movq	-96(%rbp), %rsi
	leaq	-160(%rbp), %rcx
	leaq	-180(%rbp), %rdx
	leaq	-176(%rbp), %rax
	leaq	-172(%rbp), %rbx
	movq	%rbx, 32(%rsp)
	movq	-128(%rbp), %rbx
	movq	%rbx, 24(%rsp)
	leaq	-168(%rbp), %rbx
	movq	%rbx, 16(%rsp)
	leaq	-172(%rbp), %rbx
	movq	%rbx, 8(%rsp)
	movq	-112(%rbp), %rbx
	movq	%rbx, (%rsp)
	movq	%rdi, %r9
	movq	%rsi, %r8
	movq	%rax, %rsi
	movl	$.LC29, %edi
	call	dgemv_
.LBE17:
	.loc 4 1025 0
	movl	-48(%rbp), %eax
	imull	-44(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -96(%rbp)
	.loc 4 1026 0
	cmpl	$0, -20(%rbp)
	jne	.L506
	movl	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -128(%rbp)
.L506:
	.loc 4 1012 0
	addl	$1, -76(%rbp)
.L500:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L507
	.loc 4 1028 0
	leaq	-144(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L508
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30632, %edx
	movl	$1028, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L493
.L508:
	.loc 4 1029 0
	leaq	-152(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L509
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30632, %edx
	movl	$1029, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L493
.L509:
	.loc 4 1030 0
	movq	-136(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm1
	vcvtsi2sd	-44(%rbp), %xmm0, %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	movl	-56(%rbp), %eax
	imull	-24(%rbp), %eax
	vcvtsi2sd	%eax, %xmm1, %xmm1
	vsubsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L510
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30632, %edx
	movl	$1030, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L493
.L510:
	.loc 4 1031 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L511
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L511
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L511:
	movl	$0, %eax
.L493:
	.loc 4 1032 0
	addq	$264, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE77:
	.size	MatMult_SeqBAIJ_N, .-MatMult_SeqBAIJ_N
.globl MatMultAdd_SeqBAIJ_1
	.type	MatMultAdd_SeqBAIJ_1, @function
MatMultAdd_SeqBAIJ_1:
.LFB78:
	.loc 4 1037 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$192, %rsp
	movq	%rdi, -152(%rbp)
	movq	%rsi, -160(%rbp)
	movq	%rdx, -168(%rbp)
	movq	%rcx, -176(%rbp)
	.loc 4 1038 0
	movq	-152(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 1043 0
	movq	-120(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -92(%rbp)
	movq	$0, -80(%rbp)
	movl	$0, -68(%rbp)
	.loc 4 1045 0
	movq	-120(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -44(%rbp)
	.loc 4 1047 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L545
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L546
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.30821, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1047, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L544
.L545:
	nop
	jmp	.L544
.L546:
	nop
.L544:
	.loc 4 1048 0
	leaq	-128(%rbp), %rdx
	movq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L518
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1048, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L518:
	.loc 4 1049 0
	leaq	-136(%rbp), %rdx
	movq	-168(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L520
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1049, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L520:
	.loc 4 1050 0
	movq	-176(%rbp), %rax
	cmpq	-168(%rbp), %rax
	je	.L521
	.loc 4 1051 0
	leaq	-144(%rbp), %rdx
	movq	-176(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L522
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1051, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L521:
	.loc 4 1053 0
	movq	-136(%rbp), %rax
	movq	%rax, -144(%rbp)
.L522:
	.loc 4 1056 0
	movq	-120(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1057 0
	movq	-120(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1058 0
	cmpl	$0, -44(%rbp)
	je	.L523
	.loc 4 1059 0
	movq	-176(%rbp), %rax
	cmpq	-168(%rbp), %rax
	je	.L524
	.loc 4 1060 0
	movl	-92(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-136(%rbp), %rcx
	movq	-144(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L524
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1060, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L524:
	.loc 4 1062 0
	movq	-120(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -92(%rbp)
	.loc 4 1063 0
	movq	-120(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -56(%rbp)
	.loc 4 1064 0
	movq	-120(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -80(%rbp)
	jmp	.L525
.L523:
	.loc 4 1066 0
	movq	-120(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -56(%rbp)
.L525:
	.loc 4 1069 0
	movl	$0, -88(%rbp)
	jmp	.L526
.L537:
	.loc 4 1070 0
	movq	-56(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-56(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -84(%rbp)
	.loc 4 1071 0
	addq	$4, -56(%rbp)
	.loc 4 1072 0
	cmpl	$0, -44(%rbp)
	jne	.L527
	.loc 4 1073 0
	cmpl	$0, -84(%rbp)
	setg	%al
	movzbl	%al, %eax
	addl	%eax, -68(%rbp)
	.loc 4 1074 0
	movq	-136(%rbp), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	jmp	.L528
.L527:
	.loc 4 1076 0
	movq	-136(%rbp), %rdx
	movl	-88(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-80(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
.L528:
.LBB18:
	.loc 4 1078 0
	movq	-64(%rbp), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -40(%rbp)
	movl	-84(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-64(%rbp), %rax
	movq	%rax, -32(%rbp)
	jmp	.L529
.L530:
	addq	$64, -40(%rbp)
.L529:
	movq	-40(%rbp), %rax
	cmpq	-32(%rbp), %rax
	jb	.L530
.LBE18:
.LBB19:
	.loc 4 1079 0
	movq	-104(%rbp), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	%rax, -24(%rbp)
	movl	-84(%rbp), %eax
	cltq
	salq	$4, %rax
	addq	-104(%rbp), %rax
	movq	%rax, -16(%rbp)
	jmp	.L531
.L532:
	addq	$64, -24(%rbp)
.L531:
	movq	-24(%rbp), %rax
	cmpq	-16(%rbp), %rax
	jb	.L532
.LBE19:
.LBB20:
	.loc 4 1080 0
	movl	$0, -4(%rbp)
	jmp	.L533
.L534:
	movl	-4(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-104(%rbp), %rax
	vmovsd	(%rax), %xmm1
	movq	-128(%rbp), %rdx
	movl	-4(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	-112(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -112(%rbp)
	addl	$1, -4(%rbp)
.L533:
	movl	-4(%rbp), %eax
	cmpl	-84(%rbp), %eax
	jl	.L534
.LBE20:
	.loc 4 1081 0
	movl	-84(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -104(%rbp)
	.loc 4 1082 0
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	%rax, -64(%rbp)
	.loc 4 1083 0
	cmpl	$0, -44(%rbp)
	je	.L535
	.loc 4 1084 0
	movq	-144(%rbp), %rdx
	movl	-88(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-80(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-112(%rbp), %rax
	movq	%rax, (%rdx)
	jmp	.L536
.L535:
	.loc 4 1086 0
	movq	-144(%rbp), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rdx
	movq	-112(%rbp), %rax
	movq	%rax, (%rdx)
.L536:
	.loc 4 1069 0
	addl	$1, -88(%rbp)
.L526:
	movl	-88(%rbp), %eax
	cmpl	-92(%rbp), %eax
	jl	.L537
	.loc 4 1089 0
	leaq	-128(%rbp), %rdx
	movq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L538
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1089, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L538:
	.loc 4 1090 0
	leaq	-136(%rbp), %rdx
	movq	-168(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L539
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1090, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L539:
	.loc 4 1091 0
	movq	-176(%rbp), %rax
	cmpq	-168(%rbp), %rax
	je	.L540
	.loc 4 1092 0
	leaq	-144(%rbp), %rdx
	movq	-176(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L540
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1092, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L540:
	.loc 4 1094 0
	movq	-120(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm0
	vcvtsi2sd	-68(%rbp), %xmm1, %xmm1
	vsubsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -96(%rbp)
	cmpl	$0, -96(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L541
	movl	-96(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.30821, %edx
	movl	$1094, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L519
.L541:
	.loc 4 1095 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L542
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L542
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L542:
	movl	$0, %eax
.L519:
	.loc 4 1096 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE78:
	.size	MatMultAdd_SeqBAIJ_1, .-MatMultAdd_SeqBAIJ_1
.globl MatMultAdd_SeqBAIJ_2
	.type	MatMultAdd_SeqBAIJ_2, @function
MatMultAdd_SeqBAIJ_2:
.LFB79:
	.loc 4 1101 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$240, %rsp
	movq	%rdi, -200(%rbp)
	movq	%rsi, -208(%rbp)
	movq	%rdx, -216(%rbp)
	movq	%rcx, -224(%rbp)
	.loc 4 1102 0
	movq	-200(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -160(%rbp)
	.loc 4 1103 0
	movq	$0, -152(%rbp)
	movq	$0, -144(%rbp)
	.loc 4 1107 0
	movq	-160(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	movq	$0, -48(%rbp)
	.loc 4 1108 0
	movq	-160(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 1110 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L574
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L575
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.31048, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1110, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L573
.L574:
	nop
	jmp	.L573
.L575:
	nop
.L573:
	.loc 4 1111 0
	leaq	-168(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L549
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1111, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L549:
	.loc 4 1112 0
	leaq	-176(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L551
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1112, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L551:
	.loc 4 1113 0
	movq	-224(%rbp), %rax
	cmpq	-216(%rbp), %rax
	je	.L552
	.loc 4 1114 0
	leaq	-184(%rbp), %rdx
	movq	-224(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L553
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1114, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L552:
	.loc 4 1116 0
	movq	-176(%rbp), %rax
	movq	%rax, -184(%rbp)
.L553:
	.loc 4 1119 0
	movq	-160(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1120 0
	movq	-160(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 1121 0
	cmpl	$0, -36(%rbp)
	je	.L554
	.loc 4 1122 0
	movq	-224(%rbp), %rax
	cmpq	-216(%rbp), %rax
	je	.L555
	.loc 4 1123 0
	movl	-80(%rbp), %eax
	cltq
	movq	%rax, %rdx
	salq	$4, %rdx
	movq	-176(%rbp), %rcx
	movq	-184(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L555
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1123, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L555:
	.loc 4 1125 0
	movq	-160(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1126 0
	movq	-160(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1127 0
	movq	-160(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	.loc 4 1128 0
	movq	-224(%rbp), %rax
	cmpq	-216(%rbp), %rax
	je	.L556
	.loc 4 1129 0
	movq	-160(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-176(%rbp), %rcx
	movq	-184(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L556
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1129, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L554:
	.loc 4 1132 0
	movq	-160(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1133 0
	movq	-176(%rbp), %rax
	movq	%rax, -152(%rbp)
	.loc 4 1134 0
	movq	-184(%rbp), %rax
	movq	%rax, -144(%rbp)
.L556:
	.loc 4 1137 0
	movl	$0, -76(%rbp)
	jmp	.L557
.L566:
	.loc 4 1138 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 1139 0
	cmpl	$0, -36(%rbp)
	je	.L558
	.loc 4 1140 0
	movq	-184(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$4, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -144(%rbp)
	.loc 4 1141 0
	movq	-176(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$4, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -152(%rbp)
.L558:
	.loc 4 1143 0
	movq	-152(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-152(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
.LBB21:
	.loc 4 1144 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L559
.L560:
	addq	$64, -32(%rbp)
.L559:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L560
.LBE21:
.LBB22:
	.loc 4 1145 0
	movq	-96(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$5, %rdx
	addq	%rdx, %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$6, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L561
.L562:
	addq	$64, -16(%rbp)
.L561:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L562
.LBE22:
	.loc 4 1146 0
	movl	$0, -56(%rbp)
	jmp	.L563
.L564:
	.loc 4 1147 0
	movq	-168(%rbp), %rdx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$4, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -136(%rbp)
	addq	$4, -72(%rbp)
	movq	-136(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-136(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1148 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-128(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -128(%rbp)
	.loc 4 1149 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-120(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -120(%rbp)
	.loc 4 1150 0
	addq	$32, -96(%rbp)
	.loc 4 1146 0
	addl	$1, -56(%rbp)
.L563:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L564
	.loc 4 1152 0
	movq	-144(%rbp), %rax
	movq	-128(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-144(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-120(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 1153 0
	cmpl	$0, -36(%rbp)
	jne	.L565
	.loc 4 1154 0
	addq	$16, -144(%rbp)
	addq	$16, -152(%rbp)
.L565:
	.loc 4 1137 0
	addl	$1, -76(%rbp)
.L557:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L566
	.loc 4 1157 0
	leaq	-168(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L567
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1157, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L567:
	.loc 4 1158 0
	leaq	-176(%rbp), %rdx
	movq	-216(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L568
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1158, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L568:
	.loc 4 1159 0
	movq	-224(%rbp), %rax
	cmpq	-216(%rbp), %rax
	je	.L569
	.loc 4 1160 0
	leaq	-184(%rbp), %rdx
	movq	-224(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L569
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1160, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L569:
	.loc 4 1162 0
	movq	-160(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC30(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L570
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31048, %edx
	movl	$1162, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L550
.L570:
	.loc 4 1163 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L571
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L571
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L571:
	movl	$0, %eax
.L550:
	.loc 4 1164 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE79:
	.size	MatMultAdd_SeqBAIJ_2, .-MatMultAdd_SeqBAIJ_2
.globl MatMultAdd_SeqBAIJ_3
	.type	MatMultAdd_SeqBAIJ_3, @function
MatMultAdd_SeqBAIJ_3:
.LFB80:
	.loc 4 1169 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$256, %rsp
	movq	%rdi, -216(%rbp)
	movq	%rsi, -224(%rbp)
	movq	%rdx, -232(%rbp)
	movq	%rcx, -240(%rbp)
	.loc 4 1170 0
	movq	-216(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -176(%rbp)
	.loc 4 1171 0
	movq	$0, -168(%rbp)
	movq	$0, -160(%rbp)
	.loc 4 1174 0
	movq	-176(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	movq	$0, -48(%rbp)
	.loc 4 1175 0
	movq	-176(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 1177 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L603
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L604
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.31276, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1177, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L602
.L603:
	nop
	jmp	.L602
.L604:
	nop
.L602:
	.loc 4 1178 0
	leaq	-184(%rbp), %rdx
	movq	-224(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L578
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1178, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L578:
	.loc 4 1179 0
	leaq	-192(%rbp), %rdx
	movq	-232(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L580
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1179, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L580:
	.loc 4 1180 0
	movq	-240(%rbp), %rax
	cmpq	-232(%rbp), %rax
	je	.L581
	.loc 4 1181 0
	leaq	-200(%rbp), %rdx
	movq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L582
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1181, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L581:
	.loc 4 1183 0
	movq	-192(%rbp), %rax
	movq	%rax, -200(%rbp)
.L582:
	.loc 4 1186 0
	movq	-176(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1187 0
	movq	-176(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 1188 0
	cmpl	$0, -36(%rbp)
	je	.L583
	.loc 4 1189 0
	movq	-240(%rbp), %rax
	cmpq	-232(%rbp), %rax
	je	.L584
	.loc 4 1190 0
	movl	-80(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	movq	%rax, %rdx
	movq	-192(%rbp), %rcx
	movq	-200(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L584
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1190, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L584:
	.loc 4 1192 0
	movq	-176(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1193 0
	movq	-176(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1194 0
	movq	-176(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L585
.L583:
	.loc 4 1196 0
	movq	-176(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1197 0
	movq	-192(%rbp), %rax
	movq	%rax, -168(%rbp)
	.loc 4 1198 0
	movq	-200(%rbp), %rax
	movq	%rax, -160(%rbp)
.L585:
	.loc 4 1201 0
	movl	$0, -76(%rbp)
	jmp	.L586
.L595:
	.loc 4 1202 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 1203 0
	cmpl	$0, -36(%rbp)
	je	.L587
	.loc 4 1204 0
	movq	-200(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -160(%rbp)
	.loc 4 1205 0
	movq	-192(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -168(%rbp)
.L587:
	.loc 4 1207 0
	movq	-168(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-168(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-168(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
.LBB23:
	.loc 4 1208 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L588
.L589:
	addq	$64, -32(%rbp)
.L588:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L589
.LBE23:
.LBB24:
	.loc 4 1209 0
	movq	-96(%rbp), %rcx
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L590
.L591:
	addq	$64, -16(%rbp)
.L590:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L591
.LBE24:
	.loc 4 1210 0
	movl	$0, -56(%rbp)
	jmp	.L592
.L593:
	.loc 4 1211 0
	movq	-184(%rbp), %rcx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -152(%rbp)
	addq	$4, -72(%rbp)
	movq	-152(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-152(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-152(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1212 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 1213 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 1214 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-128(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -128(%rbp)
	.loc 4 1215 0
	addq	$72, -96(%rbp)
	.loc 4 1210 0
	addl	$1, -56(%rbp)
.L592:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L593
	.loc 4 1217 0
	movq	-160(%rbp), %rax
	movq	-144(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-160(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-136(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-160(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-128(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 1218 0
	cmpl	$0, -36(%rbp)
	jne	.L594
	.loc 4 1219 0
	addq	$24, -160(%rbp)
	addq	$24, -168(%rbp)
.L594:
	.loc 4 1201 0
	addl	$1, -76(%rbp)
.L586:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L595
	.loc 4 1222 0
	leaq	-184(%rbp), %rdx
	movq	-224(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L596
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1222, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L596:
	.loc 4 1223 0
	leaq	-192(%rbp), %rdx
	movq	-232(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L597
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1223, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L597:
	.loc 4 1224 0
	movq	-240(%rbp), %rax
	cmpq	-232(%rbp), %rax
	je	.L598
	.loc 4 1225 0
	leaq	-200(%rbp), %rdx
	movq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L598
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1225, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L598:
	.loc 4 1227 0
	movq	-176(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC16(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L599
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31276, %edx
	movl	$1227, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L579
.L599:
	.loc 4 1228 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L600
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L600
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L600:
	movl	$0, %eax
.L579:
	.loc 4 1229 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE80:
	.size	MatMultAdd_SeqBAIJ_3, .-MatMultAdd_SeqBAIJ_3
.globl MatMultAdd_SeqBAIJ_4
	.type	MatMultAdd_SeqBAIJ_4, @function
MatMultAdd_SeqBAIJ_4:
.LFB81:
	.loc 4 1234 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$272, %rsp
	movq	%rdi, -232(%rbp)
	movq	%rsi, -240(%rbp)
	movq	%rdx, -248(%rbp)
	movq	%rcx, -256(%rbp)
	.loc 4 1235 0
	movq	-232(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -192(%rbp)
	.loc 4 1236 0
	movq	$0, -184(%rbp)
	movq	$0, -176(%rbp)
	.loc 4 1239 0
	movq	-192(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	movq	$0, -48(%rbp)
	.loc 4 1240 0
	movq	-192(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 1242 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L632
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L633
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.31515, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1242, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L631
.L632:
	nop
	jmp	.L631
.L633:
	nop
.L631:
	.loc 4 1243 0
	leaq	-200(%rbp), %rdx
	movq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L607
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1243, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L607:
	.loc 4 1244 0
	leaq	-208(%rbp), %rdx
	movq	-248(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L609
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1244, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L609:
	.loc 4 1245 0
	movq	-256(%rbp), %rax
	cmpq	-248(%rbp), %rax
	je	.L610
	.loc 4 1246 0
	leaq	-216(%rbp), %rdx
	movq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L611
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1246, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L610:
	.loc 4 1248 0
	movq	-208(%rbp), %rax
	movq	%rax, -216(%rbp)
.L611:
	.loc 4 1251 0
	movq	-192(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1252 0
	movq	-192(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 1253 0
	cmpl	$0, -36(%rbp)
	je	.L612
	.loc 4 1254 0
	movq	-256(%rbp), %rax
	cmpq	-248(%rbp), %rax
	je	.L613
	.loc 4 1255 0
	movl	-80(%rbp), %eax
	cltq
	movq	%rax, %rdx
	salq	$5, %rdx
	movq	-208(%rbp), %rcx
	movq	-216(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L613
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1255, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L613:
	.loc 4 1257 0
	movq	-192(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1258 0
	movq	-192(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1259 0
	movq	-192(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L614
.L612:
	.loc 4 1261 0
	movq	-192(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1262 0
	movq	-208(%rbp), %rax
	movq	%rax, -184(%rbp)
	.loc 4 1263 0
	movq	-216(%rbp), %rax
	movq	%rax, -176(%rbp)
.L614:
	.loc 4 1266 0
	movl	$0, -76(%rbp)
	jmp	.L615
.L624:
	.loc 4 1267 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 1268 0
	cmpl	$0, -36(%rbp)
	je	.L616
	.loc 4 1269 0
	movq	-216(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$5, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -176(%rbp)
	.loc 4 1270 0
	movq	-208(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$5, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -184(%rbp)
.L616:
	.loc 4 1272 0
	movq	-184(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-184(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-184(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-184(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
.LBB25:
	.loc 4 1273 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L617
.L618:
	addq	$64, -32(%rbp)
.L617:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L618
.LBE25:
.LBB26:
	.loc 4 1274 0
	movq	-96(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$7, %rdx
	addq	%rdx, %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$8, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L619
.L620:
	addq	$64, -16(%rbp)
.L619:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L620
.LBE26:
	.loc 4 1275 0
	movl	$0, -56(%rbp)
	jmp	.L621
.L622:
	.loc 4 1276 0
	movq	-200(%rbp), %rdx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$5, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -168(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 1277 0
	movq	-168(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-168(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-168(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-168(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1278 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 1279 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 1280 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 1281 0
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-136(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)
	.loc 4 1282 0
	subq	$-128, -96(%rbp)
	.loc 4 1275 0
	addl	$1, -56(%rbp)
.L621:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L622
	.loc 4 1284 0
	movq	-176(%rbp), %rax
	movq	-160(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-176(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-152(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-176(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-144(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-176(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-136(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 1285 0
	cmpl	$0, -36(%rbp)
	jne	.L623
	.loc 4 1286 0
	addq	$32, -176(%rbp)
	addq	$32, -184(%rbp)
.L623:
	.loc 4 1266 0
	addl	$1, -76(%rbp)
.L615:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L624
	.loc 4 1289 0
	leaq	-200(%rbp), %rdx
	movq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L625
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1289, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L625:
	.loc 4 1290 0
	leaq	-208(%rbp), %rdx
	movq	-248(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L626
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1290, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L626:
	.loc 4 1291 0
	movq	-256(%rbp), %rax
	cmpq	-248(%rbp), %rax
	je	.L627
	.loc 4 1292 0
	leaq	-216(%rbp), %rdx
	movq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L627
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1292, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L627:
	.loc 4 1294 0
	movq	-192(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC18(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L628
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31515, %edx
	movl	$1294, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L608
.L628:
	.loc 4 1295 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L629
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L629
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L629:
	movl	$0, %eax
.L608:
	.loc 4 1296 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE81:
	.size	MatMultAdd_SeqBAIJ_4, .-MatMultAdd_SeqBAIJ_4
.globl MatMultAdd_SeqBAIJ_5
	.type	MatMultAdd_SeqBAIJ_5, @function
MatMultAdd_SeqBAIJ_5:
.LFB82:
	.loc 4 1301 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$288, %rsp
	movq	%rdi, -248(%rbp)
	movq	%rsi, -256(%rbp)
	movq	%rdx, -264(%rbp)
	movq	%rcx, -272(%rbp)
	.loc 4 1302 0
	movq	-248(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -208(%rbp)
	.loc 4 1303 0
	movq	$0, -200(%rbp)
	movq	$0, -192(%rbp)
	.loc 4 1307 0
	movq	-208(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	movq	$0, -48(%rbp)
	.loc 4 1308 0
	movq	-208(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 1310 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L661
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L662
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.31786, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1310, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L660
.L661:
	nop
	jmp	.L660
.L662:
	nop
.L660:
	.loc 4 1311 0
	leaq	-216(%rbp), %rdx
	movq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L636
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1311, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L636:
	.loc 4 1312 0
	leaq	-224(%rbp), %rdx
	movq	-264(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L638
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1312, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L638:
	.loc 4 1313 0
	movq	-272(%rbp), %rax
	cmpq	-264(%rbp), %rax
	je	.L639
	.loc 4 1314 0
	leaq	-232(%rbp), %rdx
	movq	-272(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L640
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1314, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L639:
	.loc 4 1316 0
	movq	-224(%rbp), %rax
	movq	%rax, -232(%rbp)
.L640:
	.loc 4 1319 0
	movq	-208(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1320 0
	movq	-208(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 1321 0
	cmpl	$0, -36(%rbp)
	je	.L641
	.loc 4 1322 0
	movq	-272(%rbp), %rax
	cmpq	-264(%rbp), %rax
	je	.L642
	.loc 4 1323 0
	movl	-80(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	movq	%rax, %rdx
	movq	-224(%rbp), %rcx
	movq	-232(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L642
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1323, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L642:
	.loc 4 1325 0
	movq	-208(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1326 0
	movq	-208(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1327 0
	movq	-208(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L643
.L641:
	.loc 4 1329 0
	movq	-208(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1330 0
	movq	-224(%rbp), %rax
	movq	%rax, -200(%rbp)
	.loc 4 1331 0
	movq	-232(%rbp), %rax
	movq	%rax, -192(%rbp)
.L643:
	.loc 4 1334 0
	movl	$0, -76(%rbp)
	jmp	.L644
.L653:
	.loc 4 1335 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 1336 0
	cmpl	$0, -36(%rbp)
	je	.L645
	.loc 4 1337 0
	movq	-232(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -192(%rbp)
	.loc 4 1338 0
	movq	-224(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -200(%rbp)
.L645:
	.loc 4 1340 0
	movq	-200(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	movq	-200(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-200(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-200(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-200(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
.LBB27:
	.loc 4 1341 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L646
.L647:
	addq	$64, -32(%rbp)
.L646:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L647
.LBE27:
.LBB28:
	.loc 4 1342 0
	movq	-96(%rbp), %rcx
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	leaq	0(,%rax,4), %rdx
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	leaq	0(,%rax,4), %rdx
	addq	%rdx, %rax
	salq	$4, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L648
.L649:
	addq	$64, -16(%rbp)
.L648:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L649
.LBE28:
	.loc 4 1343 0
	movl	$0, -56(%rbp)
	jmp	.L650
.L651:
	.loc 4 1344 0
	movq	-216(%rbp), %rcx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -184(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 1345 0
	movq	-184(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-184(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-184(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-184(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-184(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1346 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 1347 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 1348 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 1349 0
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 1350 0
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-144(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)
	.loc 4 1351 0
	addq	$200, -96(%rbp)
	.loc 4 1343 0
	addl	$1, -56(%rbp)
.L650:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L651
	.loc 4 1353 0
	movq	-192(%rbp), %rax
	movq	-176(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-192(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-168(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-192(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-160(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-192(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-152(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-192(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-144(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 1354 0
	cmpl	$0, -36(%rbp)
	jne	.L652
	.loc 4 1355 0
	addq	$40, -192(%rbp)
	addq	$40, -200(%rbp)
.L652:
	.loc 4 1334 0
	addl	$1, -76(%rbp)
.L644:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L653
	.loc 4 1358 0
	leaq	-216(%rbp), %rdx
	movq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L654
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1358, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L654:
	.loc 4 1359 0
	leaq	-224(%rbp), %rdx
	movq	-264(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L655
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1359, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L655:
	.loc 4 1360 0
	movq	-272(%rbp), %rax
	cmpq	-264(%rbp), %rax
	je	.L656
	.loc 4 1361 0
	leaq	-232(%rbp), %rdx
	movq	-272(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L656
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1361, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L656:
	.loc 4 1363 0
	movq	-208(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC20(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L657
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.31786, %edx
	movl	$1363, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L637
.L657:
	.loc 4 1364 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L658
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L658
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L658:
	movl	$0, %eax
.L637:
	.loc 4 1365 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE82:
	.size	MatMultAdd_SeqBAIJ_5, .-MatMultAdd_SeqBAIJ_5
.globl MatMultAdd_SeqBAIJ_6
	.type	MatMultAdd_SeqBAIJ_6, @function
MatMultAdd_SeqBAIJ_6:
.LFB83:
	.loc 4 1369 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$304, %rsp
	movq	%rdi, -264(%rbp)
	movq	%rsi, -272(%rbp)
	movq	%rdx, -280(%rbp)
	movq	%rcx, -288(%rbp)
	.loc 4 1370 0
	movq	-264(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -224(%rbp)
	.loc 4 1371 0
	movq	$0, -216(%rbp)
	movq	$0, -208(%rbp)
	.loc 4 1375 0
	movq	-224(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	movq	$0, -48(%rbp)
	.loc 4 1376 0
	movq	-224(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 1378 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L690
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L691
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.32097, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1378, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L689
.L690:
	nop
	jmp	.L689
.L691:
	nop
.L689:
	.loc 4 1379 0
	leaq	-232(%rbp), %rdx
	movq	-272(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L665
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1379, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L665:
	.loc 4 1380 0
	leaq	-240(%rbp), %rdx
	movq	-280(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L667
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1380, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L667:
	.loc 4 1381 0
	movq	-288(%rbp), %rax
	cmpq	-280(%rbp), %rax
	je	.L668
	.loc 4 1382 0
	leaq	-248(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L669
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1382, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L668:
	.loc 4 1384 0
	movq	-240(%rbp), %rax
	movq	%rax, -248(%rbp)
.L669:
	.loc 4 1387 0
	movq	-224(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1388 0
	movq	-224(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 1389 0
	cmpl	$0, -36(%rbp)
	je	.L670
	.loc 4 1390 0
	movq	-288(%rbp), %rax
	cmpq	-280(%rbp), %rax
	je	.L671
	.loc 4 1391 0
	movl	-80(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	movq	%rax, %rdx
	movq	-240(%rbp), %rcx
	movq	-248(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L671
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1391, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L671:
	.loc 4 1393 0
	movq	-224(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1394 0
	movq	-224(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1395 0
	movq	-224(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -48(%rbp)
	jmp	.L672
.L670:
	.loc 4 1397 0
	movq	-224(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1398 0
	movq	-240(%rbp), %rax
	movq	%rax, -216(%rbp)
	.loc 4 1399 0
	movq	-248(%rbp), %rax
	movq	%rax, -208(%rbp)
.L672:
	.loc 4 1402 0
	movl	$0, -76(%rbp)
	jmp	.L673
.L682:
	.loc 4 1403 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -52(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 1404 0
	cmpl	$0, -36(%rbp)
	je	.L674
	.loc 4 1405 0
	movq	-248(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -208(%rbp)
	.loc 4 1406 0
	movq	-240(%rbp), %rcx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -216(%rbp)
.L674:
	.loc 4 1408 0
	movq	-216(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -192(%rbp)
	movq	-216(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -184(%rbp)
	movq	-216(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	movq	-216(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-216(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-216(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
.LBB29:
	.loc 4 1409 0
	movq	-72(%rbp), %rax
	movl	-52(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -32(%rbp)
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L675
.L676:
	addq	$64, -32(%rbp)
.L675:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L676
.LBE29:
.LBB30:
	.loc 4 1410 0
	movq	-96(%rbp), %rcx
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$5, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -16(%rbp)
	movl	-52(%rbp), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$3, %rax
	addq	%rdx, %rax
	salq	$6, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -8(%rbp)
	jmp	.L677
.L678:
	addq	$64, -16(%rbp)
.L677:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L678
.LBE30:
	.loc 4 1411 0
	movl	$0, -56(%rbp)
	jmp	.L679
.L680:
	.loc 4 1412 0
	movq	-232(%rbp), %rcx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -200(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 1413 0
	movq	-200(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-200(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-200(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-200(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	movq	-200(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -112(%rbp)
	movq	-200(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1414 0
	movq	-96(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 1415 0
	movq	-96(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 1416 0
	movq	-96(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 1417 0
	movq	-96(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-168(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -168(%rbp)
	.loc 4 1418 0
	movq	-96(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-160(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -160(%rbp)
	.loc 4 1419 0
	movq	-96(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm1
	movq	-96(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-112(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-96(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-104(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-152(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -152(%rbp)
	.loc 4 1420 0
	addq	$288, -96(%rbp)
	.loc 4 1411 0
	addl	$1, -56(%rbp)
.L679:
	movl	-56(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jl	.L680
	.loc 4 1422 0
	movq	-208(%rbp), %rax
	movq	-192(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-208(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-168(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-160(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-208(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-152(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 1423 0
	cmpl	$0, -36(%rbp)
	jne	.L681
	.loc 4 1424 0
	addq	$48, -208(%rbp)
	addq	$48, -216(%rbp)
.L681:
	.loc 4 1402 0
	addl	$1, -76(%rbp)
.L673:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L682
	.loc 4 1427 0
	leaq	-232(%rbp), %rdx
	movq	-272(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L683
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1427, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L683:
	.loc 4 1428 0
	leaq	-240(%rbp), %rdx
	movq	-280(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L684
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1428, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L684:
	.loc 4 1429 0
	movq	-288(%rbp), %rax
	cmpq	-280(%rbp), %rax
	je	.L685
	.loc 4 1430 0
	leaq	-248(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L685
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1430, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L685:
	.loc 4 1432 0
	movq	-224(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC22(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L686
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32097, %edx
	movl	$1432, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L666
.L686:
	.loc 4 1433 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L687
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L687
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L687:
	movl	$0, %eax
.L666:
	.loc 4 1434 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE83:
	.size	MatMultAdd_SeqBAIJ_6, .-MatMultAdd_SeqBAIJ_6
.globl MatMultAdd_SeqBAIJ_7
	.type	MatMultAdd_SeqBAIJ_7, @function
MatMultAdd_SeqBAIJ_7:
.LFB84:
	.loc 4 1439 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$328, %rsp
	movq	%rdi, -296(%rbp)
	movq	%rsi, -304(%rbp)
	movq	%rdx, -312(%rbp)
	movq	%rcx, -320(%rbp)
	.loc 4 1440 0
	movq	-296(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -256(%rbp)
	.loc 4 1441 0
	movq	$0, -248(%rbp)
	movq	$0, -240(%rbp)
	.loc 4 1445 0
	movq	-256(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -96(%rbp)
	movq	$0, -64(%rbp)
	.loc 4 1446 0
	movq	-256(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -52(%rbp)
	.loc 4 1448 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L719
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L720
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.32456, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1448, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L718
.L719:
	nop
	jmp	.L718
.L720:
	nop
.L718:
	.loc 4 1449 0
	leaq	-264(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L694
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1449, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L694:
	.loc 4 1450 0
	leaq	-272(%rbp), %rdx
	movq	-312(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L696
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1450, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L696:
	.loc 4 1451 0
	movq	-320(%rbp), %rax
	cmpq	-312(%rbp), %rax
	je	.L697
	.loc 4 1452 0
	leaq	-280(%rbp), %rdx
	movq	-320(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L698
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1452, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L697:
	.loc 4 1454 0
	movq	-272(%rbp), %rax
	movq	%rax, -280(%rbp)
.L698:
	.loc 4 1457 0
	movq	-256(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -88(%rbp)
	.loc 4 1458 0
	movq	-256(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 1459 0
	cmpl	$0, -52(%rbp)
	je	.L699
	.loc 4 1460 0
	movq	-320(%rbp), %rax
	cmpq	-312(%rbp), %rax
	je	.L700
	.loc 4 1461 0
	movl	-96(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	0(,%rax,8), %rdx
	subq	%rax, %rdx
	movq	-272(%rbp), %rcx
	movq	-280(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	PetscMemcpy
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L700
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1461, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L700:
	.loc 4 1463 0
	movq	-256(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -96(%rbp)
	.loc 4 1464 0
	movq	-256(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 1465 0
	movq	-256(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -64(%rbp)
	jmp	.L701
.L699:
	.loc 4 1467 0
	movq	-256(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 1468 0
	movq	-272(%rbp), %rax
	movq	%rax, -248(%rbp)
	.loc 4 1469 0
	movq	-280(%rbp), %rax
	movq	%rax, -240(%rbp)
.L701:
	.loc 4 1472 0
	movl	$0, -92(%rbp)
	jmp	.L702
.L711:
	.loc 4 1473 0
	movq	-80(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-80(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -68(%rbp)
	addq	$4, -80(%rbp)
	.loc 4 1474 0
	cmpl	$0, -52(%rbp)
	je	.L703
	.loc 4 1475 0
	movq	-280(%rbp), %rdx
	movl	-92(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	0(,%rax,8), %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -240(%rbp)
	.loc 4 1476 0
	movq	-272(%rbp), %rdx
	movl	-92(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	0(,%rax,8), %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -248(%rbp)
.L703:
	.loc 4 1478 0
	movq	-248(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -224(%rbp)
	movq	-248(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -216(%rbp)
	movq	-248(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -208(%rbp)
	movq	-248(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -200(%rbp)
	movq	-248(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -192(%rbp)
	movq	-248(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -184(%rbp)
	movq	-248(%rbp), %rax
	addq	$48, %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
.LBB31:
	.loc 4 1479 0
	movq	-88(%rbp), %rax
	movl	-68(%rbp), %edx
	movslq	%edx, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movq	%rax, -48(%rbp)
	movl	-68(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-88(%rbp), %rax
	movq	%rax, -40(%rbp)
	jmp	.L704
.L705:
	addq	$64, -48(%rbp)
.L704:
	movq	-48(%rbp), %rax
	cmpq	-40(%rbp), %rax
	jb	.L705
.LBE31:
.LBB32:
	.loc 4 1480 0
	movq	-112(%rbp), %rdx
	movl	-68(%rbp), %eax
	cltq
	imulq	$392, %rax, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -32(%rbp)
	movl	-68(%rbp), %eax
	cltq
	imulq	$784, %rax, %rax
	addq	-112(%rbp), %rax
	movq	%rax, -24(%rbp)
	jmp	.L706
.L707:
	addq	$64, -32(%rbp)
.L706:
	movq	-32(%rbp), %rax
	cmpq	-24(%rbp), %rax
	jb	.L707
.LBE32:
	.loc 4 1481 0
	movl	$0, -72(%rbp)
	jmp	.L708
.L709:
	.loc 4 1482 0
	movq	-264(%rbp), %rdx
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	0(,%rax,8), %rcx
	movq	%rcx, %rbx
	subq	%rax, %rbx
	movq	%rbx, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -232(%rbp)
	addq	$4, -88(%rbp)
	.loc 4 1483 0
	movq	-232(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-232(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-232(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-232(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-232(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-232(%rbp), %rax
	addq	$40, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	movq	-232(%rbp), %rax
	addq	$48, %rax
	movq	(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 1484 0
	movq	-112(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$224, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$280, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$336, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -224(%rbp)
	.loc 4 1485 0
	movq	-112(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$232, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$288, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$344, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-216(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -216(%rbp)
	.loc 4 1486 0
	movq	-112(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$240, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$296, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$352, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-208(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -208(%rbp)
	.loc 4 1487 0
	movq	-112(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$248, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$304, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$360, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-200(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -200(%rbp)
	.loc 4 1488 0
	movq	-112(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$200, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$256, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$312, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$368, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-192(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -192(%rbp)
	.loc 4 1489 0
	movq	-112(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$208, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$264, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$320, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$376, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-184(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -184(%rbp)
	.loc 4 1490 0
	movq	-112(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm1
	movq	-112(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$216, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$272, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$328, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-112(%rbp), %rax
	addq	$384, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-120(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	-176(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -176(%rbp)
	.loc 4 1491 0
	addq	$392, -112(%rbp)
	.loc 4 1481 0
	addl	$1, -72(%rbp)
.L708:
	movl	-72(%rbp), %eax
	cmpl	-68(%rbp), %eax
	jl	.L709
	.loc 4 1493 0
	movq	-240(%rbp), %rax
	movq	-224(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-240(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-216(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	16(%rax), %rdx
	movq	-208(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	24(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	32(%rax), %rdx
	movq	-192(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	40(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-240(%rbp), %rax
	leaq	48(%rax), %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
	.loc 4 1494 0
	cmpl	$0, -52(%rbp)
	jne	.L710
	.loc 4 1495 0
	addq	$56, -240(%rbp)
	addq	$56, -248(%rbp)
.L710:
	.loc 4 1472 0
	addl	$1, -92(%rbp)
.L702:
	movl	-92(%rbp), %eax
	cmpl	-96(%rbp), %eax
	jl	.L711
	.loc 4 1498 0
	leaq	-264(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L712
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1498, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L712:
	.loc 4 1499 0
	leaq	-272(%rbp), %rdx
	movq	-312(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L713
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1499, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L713:
	.loc 4 1500 0
	movq	-320(%rbp), %rax
	cmpq	-312(%rbp), %rax
	je	.L714
	.loc 4 1501 0
	leaq	-280(%rbp), %rdx
	movq	-320(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L714
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1501, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L714:
	.loc 4 1503 0
	movq	-256(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	.LC24(%rip), %xmm1
	vmulsd	%xmm1, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -100(%rbp)
	cmpl	$0, -100(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L715
	movl	-100(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32456, %edx
	movl	$1503, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L695
.L715:
	.loc 4 1504 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L716
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L716
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L716:
	movl	$0, %eax
.L695:
	.loc 4 1505 0
	addq	$328, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE84:
	.size	MatMultAdd_SeqBAIJ_7, .-MatMultAdd_SeqBAIJ_7
.globl MatMultAdd_SeqBAIJ_N
	.type	MatMultAdd_SeqBAIJ_N, @function
MatMultAdd_SeqBAIJ_N:
.LFB85:
	.loc 4 1510 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$248, %rsp
	movq	%rdi, -184(%rbp)
	movq	%rsi, -192(%rbp)
	movq	%rdx, -200(%rbp)
	movq	%rcx, -208(%rbp)
	.loc 4 1511 0
	movq	-184(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -136(%rbp)
	.loc 4 1512 0
	movq	$0, -128(%rbp)
	.loc 4 1515 0
	movq	-184(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -56(%rbp)
	movq	-136(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -44(%rbp)
	.loc 4 1516 0
	movq	$0, -32(%rbp)
	.loc 4 1517 0
	movq	-136(%rbp), %rax
	movl	100(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 1519 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L746
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L747
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.32859, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1519, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L745
.L746:
	nop
	jmp	.L745
.L747:
	nop
.L745:
	.loc 4 1520 0
	movq	-208(%rbp), %rdx
	movq	-200(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecCopy
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L723
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1520, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L724
.L723:
	.loc 4 1521 0
	leaq	-144(%rbp), %rdx
	movq	-192(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L725
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1521, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L724
.L725:
	.loc 4 1522 0
	leaq	-152(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L726
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1522, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L724
.L726:
	.loc 4 1524 0
	movq	-136(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1525 0
	movq	-136(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 1526 0
	cmpl	$0, -20(%rbp)
	je	.L727
	.loc 4 1527 0
	movq	-136(%rbp), %rax
	movl	104(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1528 0
	movq	-136(%rbp), %rax
	movq	112(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1529 0
	movq	-136(%rbp), %rax
	movq	120(%rax), %rax
	movq	%rax, -32(%rbp)
	jmp	.L728
.L727:
	.loc 4 1531 0
	movq	-136(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -80(%rbp)
	.loc 4 1532 0
	movq	-136(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -64(%rbp)
	.loc 4 1533 0
	movq	-152(%rbp), %rax
	movq	%rax, -128(%rbp)
.L728:
	.loc 4 1536 0
	movq	-136(%rbp), %rax
	movq	240(%rax), %rax
	testq	%rax, %rax
	jne	.L729
	.loc 4 1537 0
	movq	-184(%rbp), %rax
	movq	456(%rax), %rax
	movl	4(%rax), %edx
	movq	-184(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	cmpl	%eax, %edx
	cmovge	%edx, %eax
	movl	%eax, -36(%rbp)
	.loc 4 1538 0
	movl	-36(%rbp), %eax
	addl	$1, %eax
	cltq
	salq	$3, %rax
	testq	%rax, %rax
	je	.L730
	movq	PetscTrMalloc(%rip), %rax
	movq	-136(%rbp), %rdx
	addq	$240, %rdx
	movl	-36(%rbp), %ecx
	addl	$1, %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,8), %rbx
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1538, %esi
	movq	%rbx, %rdi
	call	*%rax
	jmp	.L731
.L730:
	movq	-136(%rbp), %rax
	movq	$0, 240(%rax)
	movl	$0, %eax
.L731:
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L729
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1538, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L724
.L729:
	.loc 4 1540 0
	movq	-136(%rbp), %rax
	movq	240(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 1541 0
	movl	$0, -76(%rbp)
	jmp	.L732
.L739:
	.loc 4 1542 0
	movq	-64(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-64(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -48(%rbp)
	addq	$4, -64(%rbp)
	.loc 4 1543 0
	movl	-48(%rbp), %eax
	imull	-56(%rbp), %eax
	movl	%eax, -40(%rbp)
	.loc 4 1544 0
	movq	-112(%rbp), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1545 0
	movl	$0, -52(%rbp)
	jmp	.L733
.L736:
	.loc 4 1546 0
	movq	-144(%rbp), %rdx
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	imull	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -120(%rbp)
	addq	$4, -72(%rbp)
	.loc 4 1547 0
	movl	$0, -36(%rbp)
	jmp	.L734
.L735:
	movl	-36(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-104(%rbp), %rax
	movl	-36(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	-120(%rbp), %rdx
	movq	(%rdx), %rdx
	movq	%rdx, (%rax)
	addl	$1, -36(%rbp)
.L734:
	movl	-36(%rbp), %eax
	cmpl	-56(%rbp), %eax
	jl	.L735
	.loc 4 1548 0
	movl	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -104(%rbp)
	.loc 4 1545 0
	addl	$1, -52(%rbp)
.L733:
	movl	-52(%rbp), %eax
	cmpl	-48(%rbp), %eax
	jl	.L736
	.loc 4 1550 0
	cmpl	$0, -20(%rbp)
	je	.L737
	movq	-152(%rbp), %rdx
	movl	-76(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-32(%rbp), %rax
	movl	(%rax), %eax
	imull	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -128(%rbp)
.L737:
.LBB33:
	.loc 4 1551 0
	movabsq	$4607182418800017408, %rax
	movq	%rax, -160(%rbp)
	movl	$1, -164(%rbp)
	movl	-56(%rbp), %eax
	movl	%eax, -168(%rbp)
	movl	-40(%rbp), %eax
	movl	%eax, -172(%rbp)
	leaq	-168(%rbp), %rdi
	movq	-96(%rbp), %rsi
	leaq	-160(%rbp), %rcx
	leaq	-172(%rbp), %rdx
	leaq	-168(%rbp), %rax
	leaq	-164(%rbp), %rbx
	movq	%rbx, 32(%rsp)
	movq	-128(%rbp), %rbx
	movq	%rbx, 24(%rsp)
	leaq	-160(%rbp), %rbx
	movq	%rbx, 16(%rsp)
	leaq	-164(%rbp), %rbx
	movq	%rbx, 8(%rsp)
	movq	-112(%rbp), %rbx
	movq	%rbx, (%rsp)
	movq	%rdi, %r9
	movq	%rsi, %r8
	movq	%rax, %rsi
	movl	$.LC29, %edi
	call	dgemv_
.LBE33:
	.loc 4 1553 0
	movl	-48(%rbp), %eax
	imull	-44(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -96(%rbp)
	.loc 4 1554 0
	cmpl	$0, -20(%rbp)
	jne	.L738
	.loc 4 1555 0
	movl	-56(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -128(%rbp)
.L738:
	.loc 4 1541 0
	addl	$1, -76(%rbp)
.L732:
	movl	-76(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L739
	.loc 4 1558 0
	leaq	-144(%rbp), %rdx
	movq	-192(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L740
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1558, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L724
.L740:
	.loc 4 1559 0
	leaq	-152(%rbp), %rdx
	movq	-208(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L741
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1559, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L724
.L741:
	.loc 4 1560 0
	movq	-136(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm1
	vcvtsi2sd	-44(%rbp), %xmm0, %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L742
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.32859, %edx
	movl	$1560, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L724
.L742:
	.loc 4 1561 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L743
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L743
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L743:
	movl	$0, %eax
.L724:
	.loc 4 1562 0
	addq	$248, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE85:
	.size	MatMultAdd_SeqBAIJ_N, .-MatMultAdd_SeqBAIJ_N
.globl MatMultHermitianTranspose_SeqBAIJ
	.type	MatMultHermitianTranspose_SeqBAIJ, @function
MatMultHermitianTranspose_SeqBAIJ:
.LFB86:
	.loc 4 1567 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$72, %rsp
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	.loc 4 1568 0
	movl	$0, %eax
	movq	%rax, -32(%rbp)
	.loc 4 1571 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L756
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L757
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.33035, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1571, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L755
.L756:
	nop
	jmp	.L755
.L757:
	nop
.L755:
	.loc 4 1572 0
	vmovsd	-32(%rbp), %xmm0
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	call	VecSet
	movl	%eax, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L750
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33035, %edx
	movl	$1572, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L751
.L750:
	.loc 4 1573 0
	movq	-56(%rbp), %rcx
	movq	-56(%rbp), %rdx
	movq	-48(%rbp), %rbx
	movq	-40(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	MatMultHermitianTransposeAdd_SeqBAIJ
	movl	%eax, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L752
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33035, %edx
	movl	$1573, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L751
.L752:
	.loc 4 1574 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L753
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L753
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L753:
	movl	$0, %eax
.L751:
	.loc 4 1575 0
	addq	$72, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE86:
	.size	MatMultHermitianTranspose_SeqBAIJ, .-MatMultHermitianTranspose_SeqBAIJ
.globl MatMultTranspose_SeqBAIJ
	.type	MatMultTranspose_SeqBAIJ, @function
MatMultTranspose_SeqBAIJ:
.LFB87:
	.loc 4 1580 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$72, %rsp
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	.loc 4 1581 0
	movl	$0, %eax
	movq	%rax, -32(%rbp)
	.loc 4 1584 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L766
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L767
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.33103, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1584, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L765
.L766:
	nop
	jmp	.L765
.L767:
	nop
.L765:
	.loc 4 1585 0
	vmovsd	-32(%rbp), %xmm0
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	call	VecSet
	movl	%eax, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L760
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33103, %edx
	movl	$1585, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L761
.L760:
	.loc 4 1586 0
	movq	-56(%rbp), %rcx
	movq	-56(%rbp), %rdx
	movq	-48(%rbp), %rbx
	movq	-40(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	MatMultTransposeAdd_SeqBAIJ
	movl	%eax, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L762
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33103, %edx
	movl	$1586, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L761
.L762:
	.loc 4 1587 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L763
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L763
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L763:
	movl	$0, %eax
.L761:
	.loc 4 1588 0
	addq	$72, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE87:
	.size	MatMultTranspose_SeqBAIJ, .-MatMultTranspose_SeqBAIJ
	.section	.rodata
	.align 8
.LC31:
	.string	"block size larger than 5 is not supported yet"
	.text
.globl MatMultHermitianTransposeAdd_SeqBAIJ
	.type	MatMultHermitianTransposeAdd_SeqBAIJ, @function
MatMultHermitianTransposeAdd_SeqBAIJ:
.LFB88:
	.loc 4 1594 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$320, %rsp
	movq	%rdi, -280(%rbp)
	movq	%rsi, -288(%rbp)
	movq	%rdx, -296(%rbp)
	movq	%rcx, -304(%rbp)
	.loc 4 1595 0
	movq	-280(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -184(%rbp)
	.loc 4 1596 0
	movq	$0, -168(%rbp)
	.loc 4 1599 0
	movq	-280(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -72(%rbp)
	movq	-184(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -60(%rbp)
	movq	$0, -48(%rbp)
	.loc 4 1600 0
	movq	-184(%rbp), %rax
	movq	96(%rax), %rdx
	movq	%rdx, -240(%rbp)
	movq	104(%rax), %rdx
	movq	%rdx, -232(%rbp)
	movq	112(%rax), %rdx
	movq	%rdx, -224(%rbp)
	movq	120(%rax), %rax
	movq	%rax, -216(%rbp)
	.loc 4 1601 0
	movl	-236(%rbp), %eax
	movl	%eax, -36(%rbp)
	.loc 4 1603 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L820
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L821
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.33195, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1603, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L819
.L820:
	nop
	jmp	.L819
.L821:
	nop
.L819:
	.loc 4 1604 0
	movq	-296(%rbp), %rax
	cmpq	-304(%rbp), %rax
	je	.L770
	movq	-304(%rbp), %rdx
	movq	-296(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecCopy
	movl	%eax, -108(%rbp)
	cmpl	$0, -108(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L770
	movl	-108(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33195, %edx
	movl	$1604, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L771
.L770:
	.loc 4 1605 0
	leaq	-192(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -108(%rbp)
	cmpl	$0, -108(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L772
	movl	-108(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33195, %edx
	movl	$1605, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L771
.L772:
	.loc 4 1606 0
	leaq	-200(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -108(%rbp)
	cmpl	$0, -108(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L773
	movl	-108(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33195, %edx
	movl	$1606, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L771
.L773:
	.loc 4 1608 0
	movq	-184(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 1609 0
	movq	-184(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -120(%rbp)
	.loc 4 1610 0
	cmpl	$0, -36(%rbp)
	je	.L774
	.loc 4 1611 0
	movl	-232(%rbp), %eax
	movl	%eax, -104(%rbp)
	.loc 4 1612 0
	movq	-224(%rbp), %rax
	movq	%rax, -88(%rbp)
	.loc 4 1613 0
	movq	-216(%rbp), %rax
	movq	%rax, -48(%rbp)
	jmp	.L775
.L774:
	.loc 4 1615 0
	movq	-184(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -104(%rbp)
	.loc 4 1616 0
	movq	-184(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -88(%rbp)
	.loc 4 1617 0
	movq	-192(%rbp), %rax
	movq	%rax, -168(%rbp)
.L775:
	.loc 4 1620 0
	cmpl	$5, -72(%rbp)
	ja	.L776
	mov	-72(%rbp), %eax
	movq	.L782(,%rax,8), %rax
	jmp	*%rax
	.section	.rodata
	.align 8
	.align 4
.L782:
	.quad	.L776
	.quad	.L777
	.quad	.L778
	.quad	.L779
	.quad	.L780
	.quad	.L781
	.text
.L777:
	.loc 4 1622 0
	movl	$0, -100(%rbp)
	jmp	.L783
.L788:
	.loc 4 1623 0
	cmpl	$0, -36(%rbp)
	je	.L784
	movq	-192(%rbp), %rdx
	movl	-100(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -168(%rbp)
.L784:
	.loc 4 1624 0
	movq	-168(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	.loc 4 1625 0
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 1626 0
	movq	-88(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -64(%rbp)
	addq	$4, -88(%rbp)
	.loc 4 1627 0
	movl	$0, -68(%rbp)
	jmp	.L785
.L786:
	.loc 4 1628 0
	movl	-68(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -76(%rbp)
	.loc 4 1629 0
	movq	-200(%rbp), %rax
	movl	-76(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rdx
	movq	-200(%rbp), %rax
	movl	-76(%rbp), %ecx
	movslq	%ecx, %rcx
	salq	$3, %rcx
	addq	%rcx, %rax
	vmovsd	(%rax), %xmm1
	movq	-120(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, (%rdx)
	.loc 4 1630 0
	addq	$8, -120(%rbp)
	.loc 4 1627 0
	addl	$1, -68(%rbp)
.L785:
	movl	-68(%rbp), %eax
	cmpl	-64(%rbp), %eax
	jl	.L786
	.loc 4 1632 0
	cmpl	$0, -36(%rbp)
	jne	.L787
	addq	$8, -168(%rbp)
.L787:
	.loc 4 1622 0
	addl	$1, -100(%rbp)
.L783:
	movl	-100(%rbp), %eax
	cmpl	-104(%rbp), %eax
	jl	.L788
	.loc 4 1634 0
	jmp	.L789
.L778:
	.loc 4 1636 0
	movl	$0, -100(%rbp)
	jmp	.L790
.L795:
	.loc 4 1637 0
	cmpl	$0, -36(%rbp)
	je	.L791
	movq	-192(%rbp), %rdx
	movl	-100(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$4, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -168(%rbp)
.L791:
	.loc 4 1638 0
	movq	-168(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-168(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	.loc 4 1639 0
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 1640 0
	movq	-88(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -64(%rbp)
	addq	$4, -88(%rbp)
	.loc 4 1641 0
	movl	$0, -68(%rbp)
	jmp	.L792
.L793:
	.loc 4 1642 0
	movl	-68(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	addl	%eax, %eax
	movl	%eax, -76(%rbp)
	.loc 4 1643 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1644 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1645 0
	addq	$32, -120(%rbp)
	.loc 4 1641 0
	addl	$1, -68(%rbp)
.L792:
	movl	-68(%rbp), %eax
	cmpl	-64(%rbp), %eax
	jl	.L793
	.loc 4 1647 0
	cmpl	$0, -36(%rbp)
	jne	.L794
	addq	$16, -168(%rbp)
.L794:
	.loc 4 1636 0
	addl	$1, -100(%rbp)
.L790:
	movl	-100(%rbp), %eax
	cmpl	-104(%rbp), %eax
	jl	.L795
	.loc 4 1649 0
	jmp	.L789
.L779:
	.loc 4 1651 0
	movl	$0, -100(%rbp)
	jmp	.L796
.L801:
	.loc 4 1652 0
	cmpl	$0, -36(%rbp)
	je	.L797
	movq	-192(%rbp), %rcx
	movl	-100(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -168(%rbp)
.L797:
	.loc 4 1653 0
	movq	-168(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-168(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-168(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	.loc 4 1654 0
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 1655 0
	movq	-88(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -64(%rbp)
	addq	$4, -88(%rbp)
	.loc 4 1656 0
	movl	$0, -68(%rbp)
	jmp	.L798
.L799:
	.loc 4 1657 0
	movl	-68(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %edx
	movl	%edx, %eax
	addl	%eax, %eax
	addl	%edx, %eax
	movl	%eax, -76(%rbp)
	.loc 4 1658 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1659 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1660 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1661 0
	addq	$72, -120(%rbp)
	.loc 4 1656 0
	addl	$1, -68(%rbp)
.L798:
	movl	-68(%rbp), %eax
	cmpl	-64(%rbp), %eax
	jl	.L799
	.loc 4 1663 0
	cmpl	$0, -36(%rbp)
	jne	.L800
	addq	$24, -168(%rbp)
.L800:
	.loc 4 1651 0
	addl	$1, -100(%rbp)
.L796:
	movl	-100(%rbp), %eax
	cmpl	-104(%rbp), %eax
	jl	.L801
	.loc 4 1665 0
	jmp	.L789
.L780:
	.loc 4 1667 0
	movl	$0, -100(%rbp)
	jmp	.L802
.L807:
	.loc 4 1668 0
	cmpl	$0, -36(%rbp)
	je	.L803
	movq	-192(%rbp), %rdx
	movl	-100(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$5, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -168(%rbp)
.L803:
	.loc 4 1669 0
	movq	-168(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-168(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-168(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	movq	-168(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	.loc 4 1670 0
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 1671 0
	movq	-88(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -64(%rbp)
	addq	$4, -88(%rbp)
	.loc 4 1672 0
	movl	$0, -68(%rbp)
	jmp	.L804
.L805:
	.loc 4 1673 0
	movl	-68(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %eax
	sall	$2, %eax
	movl	%eax, -76(%rbp)
	.loc 4 1674 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1675 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1676 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1677 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1678 0
	subq	$-128, -120(%rbp)
	.loc 4 1672 0
	addl	$1, -68(%rbp)
.L804:
	movl	-68(%rbp), %eax
	cmpl	-64(%rbp), %eax
	jl	.L805
	.loc 4 1680 0
	cmpl	$0, -36(%rbp)
	jne	.L806
	addq	$32, -168(%rbp)
.L806:
	.loc 4 1667 0
	addl	$1, -100(%rbp)
.L802:
	movl	-100(%rbp), %eax
	cmpl	-104(%rbp), %eax
	jl	.L807
	.loc 4 1682 0
	jmp	.L789
.L781:
	.loc 4 1684 0
	movl	$0, -100(%rbp)
	jmp	.L808
.L813:
	.loc 4 1685 0
	cmpl	$0, -36(%rbp)
	je	.L809
	movq	-192(%rbp), %rcx
	movl	-100(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -168(%rbp)
.L809:
	.loc 4 1686 0
	movq	-168(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-168(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-168(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	.loc 4 1687 0
	movq	-168(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -136(%rbp)
	movq	-168(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -128(%rbp)
	.loc 4 1688 0
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-96(%rbp), %rax
	movq	%rax, -56(%rbp)
	.loc 4 1689 0
	movq	-88(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-88(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -64(%rbp)
	addq	$4, -88(%rbp)
	.loc 4 1690 0
	movl	$0, -68(%rbp)
	jmp	.L810
.L811:
	.loc 4 1691 0
	movl	-68(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-56(%rbp), %rax
	movl	(%rax), %edx
	movl	%edx, %eax
	sall	$2, %eax
	addl	%edx, %eax
	movl	%eax, -76(%rbp)
	.loc 4 1692 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1693 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1694 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1695 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1696 0
	movq	-200(%rbp), %rdx
	movl	-76(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-200(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-120(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm1
	movq	-120(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-136(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-120(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-128(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -76(%rbp)
	.loc 4 1697 0
	addq	$200, -120(%rbp)
	.loc 4 1690 0
	addl	$1, -68(%rbp)
.L810:
	movl	-68(%rbp), %eax
	cmpl	-64(%rbp), %eax
	jl	.L811
	.loc 4 1699 0
	cmpl	$0, -36(%rbp)
	jne	.L812
	addq	$40, -168(%rbp)
.L812:
	.loc 4 1684 0
	addl	$1, -100(%rbp)
.L808:
	movl	-100(%rbp), %eax
	cmpl	-104(%rbp), %eax
	jl	.L813
	.loc 4 1701 0
	jmp	.L789
.L776:
.LBB34:
	.loc 4 1706 0
	movq	$.LC31, 8(%rsp)
	movl	$0, (%rsp)
	movl	$56, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33195, %edx
	movl	$1706, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L771
.L789:
.LBE34:
	.loc 4 1733 0
	leaq	-192(%rbp), %rdx
	movq	-288(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -108(%rbp)
	cmpl	$0, -108(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L814
	movl	-108(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33195, %edx
	movl	$1733, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L771
.L814:
	.loc 4 1734 0
	leaq	-200(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -108(%rbp)
	cmpl	$0, -108(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L815
	movl	-108(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33195, %edx
	movl	$1734, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L771
.L815:
	.loc 4 1735 0
	movq	-184(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm1
	movq	-184(%rbp), %rax
	movl	224(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -108(%rbp)
	cmpl	$0, -108(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L816
	movl	-108(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33195, %edx
	movl	$1735, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L771
.L816:
	.loc 4 1736 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L817
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L817
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L817:
	movl	$0, %eax
.L771:
	.loc 4 1737 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE88:
	.size	MatMultHermitianTransposeAdd_SeqBAIJ, .-MatMultHermitianTransposeAdd_SeqBAIJ
	.section	.rodata
.LC32:
	.string	"T"
	.text
.globl MatMultTransposeAdd_SeqBAIJ
	.type	MatMultTransposeAdd_SeqBAIJ, @function
MatMultTransposeAdd_SeqBAIJ:
.LFB89:
	.loc 4 1742 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$360, %rsp
	movq	%rdi, -296(%rbp)
	movq	%rsi, -304(%rbp)
	movq	%rdx, -312(%rbp)
	movq	%rcx, -320(%rbp)
	.loc 4 1743 0
	movq	-296(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -200(%rbp)
	.loc 4 1744 0
	movq	$0, -184(%rbp)
	.loc 4 1747 0
	movq	-296(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -88(%rbp)
	movq	-200(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -76(%rbp)
	movq	$0, -64(%rbp)
	.loc 4 1748 0
	movq	-200(%rbp), %rax
	movq	96(%rax), %rdx
	movq	%rdx, -256(%rbp)
	movq	104(%rax), %rdx
	movq	%rdx, -248(%rbp)
	movq	112(%rax), %rdx
	movq	%rdx, -240(%rbp)
	movq	120(%rax), %rax
	movq	%rax, -232(%rbp)
	.loc 4 1749 0
	movl	-252(%rbp), %eax
	movl	%eax, -52(%rbp)
	.loc 4 1751 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L886
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L887
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.33932, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1751, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L885
.L886:
	nop
	jmp	.L885
.L887:
	nop
.L885:
	.loc 4 1752 0
	movq	-312(%rbp), %rax
	cmpq	-320(%rbp), %rax
	je	.L824
	movq	-320(%rbp), %rdx
	movq	-312(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecCopy
	movl	%eax, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L824
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1752, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L824:
	.loc 4 1753 0
	leaq	-208(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L826
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1753, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L826:
	.loc 4 1754 0
	leaq	-216(%rbp), %rdx
	movq	-320(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L827
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1754, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L827:
	.loc 4 1756 0
	movq	-200(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 1757 0
	movq	-200(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -136(%rbp)
	.loc 4 1758 0
	cmpl	$0, -52(%rbp)
	je	.L828
	.loc 4 1759 0
	movl	-248(%rbp), %eax
	movl	%eax, -120(%rbp)
	.loc 4 1760 0
	movq	-240(%rbp), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1761 0
	movq	-232(%rbp), %rax
	movq	%rax, -64(%rbp)
	jmp	.L829
.L828:
	.loc 4 1763 0
	movq	-200(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -120(%rbp)
	.loc 4 1764 0
	movq	-200(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 1765 0
	movq	-208(%rbp), %rax
	movq	%rax, -184(%rbp)
.L829:
	.loc 4 1768 0
	cmpl	$5, -88(%rbp)
	ja	.L830
	mov	-88(%rbp), %eax
	movq	.L836(,%rax,8), %rax
	jmp	*%rax
	.section	.rodata
	.align 8
	.align 4
.L836:
	.quad	.L830
	.quad	.L831
	.quad	.L832
	.quad	.L833
	.quad	.L834
	.quad	.L835
	.text
.L831:
	.loc 4 1770 0
	movl	$0, -116(%rbp)
	jmp	.L837
.L842:
	.loc 4 1771 0
	cmpl	$0, -52(%rbp)
	je	.L838
	movq	-208(%rbp), %rdx
	movl	-116(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -184(%rbp)
.L838:
	.loc 4 1772 0
	movq	-184(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	.loc 4 1773 0
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-112(%rbp), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1774 0
	movq	-104(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -80(%rbp)
	addq	$4, -104(%rbp)
	.loc 4 1775 0
	movl	$0, -84(%rbp)
	jmp	.L839
.L840:
	.loc 4 1776 0
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -92(%rbp)
	.loc 4 1777 0
	movq	-216(%rbp), %rax
	movl	-92(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rdx
	movq	-216(%rbp), %rax
	movl	-92(%rbp), %ecx
	movslq	%ecx, %rcx
	salq	$3, %rcx
	addq	%rcx, %rax
	vmovsd	(%rax), %xmm1
	movq	-136(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, (%rdx)
	.loc 4 1778 0
	addq	$8, -136(%rbp)
	.loc 4 1775 0
	addl	$1, -84(%rbp)
.L839:
	movl	-84(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L840
	.loc 4 1780 0
	cmpl	$0, -52(%rbp)
	jne	.L841
	addq	$8, -184(%rbp)
.L841:
	.loc 4 1770 0
	addl	$1, -116(%rbp)
.L837:
	movl	-116(%rbp), %eax
	cmpl	-120(%rbp), %eax
	jl	.L842
	.loc 4 1782 0
	jmp	.L843
.L832:
	.loc 4 1784 0
	movl	$0, -116(%rbp)
	jmp	.L844
.L849:
	.loc 4 1785 0
	cmpl	$0, -52(%rbp)
	je	.L845
	movq	-208(%rbp), %rdx
	movl	-116(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$4, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -184(%rbp)
.L845:
	.loc 4 1786 0
	movq	-184(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	movq	-184(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	.loc 4 1787 0
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-112(%rbp), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1788 0
	movq	-104(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -80(%rbp)
	addq	$4, -104(%rbp)
	.loc 4 1789 0
	movl	$0, -84(%rbp)
	jmp	.L846
.L847:
	.loc 4 1790 0
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	addl	%eax, %eax
	movl	%eax, -92(%rbp)
	.loc 4 1791 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1792 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1793 0
	addq	$32, -136(%rbp)
	.loc 4 1789 0
	addl	$1, -84(%rbp)
.L846:
	movl	-84(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L847
	.loc 4 1795 0
	cmpl	$0, -52(%rbp)
	jne	.L848
	addq	$16, -184(%rbp)
.L848:
	.loc 4 1784 0
	addl	$1, -116(%rbp)
.L844:
	movl	-116(%rbp), %eax
	cmpl	-120(%rbp), %eax
	jl	.L849
	.loc 4 1797 0
	jmp	.L843
.L833:
	.loc 4 1799 0
	movl	$0, -116(%rbp)
	jmp	.L850
.L855:
	.loc 4 1800 0
	cmpl	$0, -52(%rbp)
	je	.L851
	movq	-208(%rbp), %rcx
	movl	-116(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -184(%rbp)
.L851:
	.loc 4 1801 0
	movq	-184(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	movq	-184(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-184(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	.loc 4 1802 0
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-112(%rbp), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1803 0
	movq	-104(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -80(%rbp)
	addq	$4, -104(%rbp)
	.loc 4 1804 0
	movl	$0, -84(%rbp)
	jmp	.L852
.L853:
	.loc 4 1805 0
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %edx
	movl	%edx, %eax
	addl	%eax, %eax
	addl	%edx, %eax
	movl	%eax, -92(%rbp)
	.loc 4 1806 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1807 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1808 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1809 0
	addq	$72, -136(%rbp)
	.loc 4 1804 0
	addl	$1, -84(%rbp)
.L852:
	movl	-84(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L853
	.loc 4 1811 0
	cmpl	$0, -52(%rbp)
	jne	.L854
	addq	$24, -184(%rbp)
.L854:
	.loc 4 1799 0
	addl	$1, -116(%rbp)
.L850:
	movl	-116(%rbp), %eax
	cmpl	-120(%rbp), %eax
	jl	.L855
	.loc 4 1813 0
	jmp	.L843
.L834:
	.loc 4 1815 0
	movl	$0, -116(%rbp)
	jmp	.L856
.L861:
	.loc 4 1816 0
	cmpl	$0, -52(%rbp)
	je	.L857
	movq	-208(%rbp), %rdx
	movl	-116(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$5, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -184(%rbp)
.L857:
	.loc 4 1817 0
	movq	-184(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	movq	-184(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-184(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	movq	-184(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	.loc 4 1818 0
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-112(%rbp), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1819 0
	movq	-104(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -80(%rbp)
	addq	$4, -104(%rbp)
	.loc 4 1820 0
	movl	$0, -84(%rbp)
	jmp	.L858
.L859:
	.loc 4 1821 0
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %eax
	sall	$2, %eax
	movl	%eax, -92(%rbp)
	.loc 4 1822 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1823 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1824 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1825 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1826 0
	subq	$-128, -136(%rbp)
	.loc 4 1820 0
	addl	$1, -84(%rbp)
.L858:
	movl	-84(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L859
	.loc 4 1828 0
	cmpl	$0, -52(%rbp)
	jne	.L860
	addq	$32, -184(%rbp)
.L860:
	.loc 4 1815 0
	addl	$1, -116(%rbp)
.L856:
	movl	-116(%rbp), %eax
	cmpl	-120(%rbp), %eax
	jl	.L861
	.loc 4 1830 0
	jmp	.L843
.L835:
	.loc 4 1832 0
	movl	$0, -116(%rbp)
	jmp	.L862
.L867:
	.loc 4 1833 0
	cmpl	$0, -52(%rbp)
	je	.L863
	movq	-208(%rbp), %rcx
	movl	-116(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	movslq	%eax, %rdx
	movq	%rdx, %rax
	salq	$2, %rax
	addq	%rdx, %rax
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	movq	%rax, -184(%rbp)
.L863:
	.loc 4 1834 0
	movq	-184(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -176(%rbp)
	movq	-184(%rbp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, -168(%rbp)
	movq	-184(%rbp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, -160(%rbp)
	.loc 4 1835 0
	movq	-184(%rbp), %rax
	addq	$24, %rax
	movq	(%rax), %rax
	movq	%rax, -152(%rbp)
	movq	-184(%rbp), %rax
	addq	$32, %rax
	movq	(%rax), %rax
	movq	%rax, -144(%rbp)
	.loc 4 1836 0
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	cltq
	salq	$2, %rax
	addq	-112(%rbp), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1837 0
	movq	-104(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -80(%rbp)
	addq	$4, -104(%rbp)
	.loc 4 1838 0
	movl	$0, -84(%rbp)
	jmp	.L864
.L865:
	.loc 4 1839 0
	movl	-84(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-72(%rbp), %rax
	movl	(%rax), %edx
	movl	%edx, %eax
	sall	$2, %eax
	addl	%edx, %eax
	movl	%eax, -92(%rbp)
	.loc 4 1840 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$8, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$16, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$24, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$32, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1841 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$40, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$48, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$56, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$64, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$72, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1842 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$80, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$88, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$96, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$104, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$112, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1843 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$120, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	subq	$-128, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$136, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$144, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$152, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1844 0
	movq	-216(%rbp), %rdx
	movl	-92(%rbp), %eax
	movslq	%eax, %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movq	-216(%rbp), %rcx
	cltq
	salq	$3, %rax
	leaq	(%rcx,%rax), %rax
	vmovsd	(%rax), %xmm2
	movq	-136(%rbp), %rax
	addq	$160, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-176(%rbp), %xmm0, %xmm1
	movq	-136(%rbp), %rax
	addq	$168, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-168(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$176, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-160(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$184, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-152(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm1
	movq	-136(%rbp), %rax
	addq	$192, %rax
	vmovsd	(%rax), %xmm0
	vmulsd	-144(%rbp), %xmm0, %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vaddsd	%xmm0, %xmm2, %xmm0
	vmovsd	%xmm0, (%rdx)
	addl	$1, -92(%rbp)
	.loc 4 1845 0
	addq	$200, -136(%rbp)
	.loc 4 1838 0
	addl	$1, -84(%rbp)
.L864:
	movl	-84(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L865
	.loc 4 1847 0
	cmpl	$0, -52(%rbp)
	jne	.L866
	addq	$40, -184(%rbp)
.L866:
	.loc 4 1832 0
	addl	$1, -116(%rbp)
.L862:
	movl	-116(%rbp), %eax
	cmpl	-120(%rbp), %eax
	jl	.L867
	.loc 4 1849 0
	jmp	.L843
.L830:
.LBB35:
	.loc 4 1854 0
	movq	-200(%rbp), %rax
	movq	240(%rax), %rax
	testq	%rax, %rax
	jne	.L868
	.loc 4 1855 0
	movq	-296(%rbp), %rax
	movq	456(%rax), %rax
	movl	4(%rax), %edx
	movq	-296(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	cmpl	%eax, %edx
	cmovge	%edx, %eax
	movl	%eax, -44(%rbp)
	.loc 4 1856 0
	movl	-44(%rbp), %eax
	addl	$1, %eax
	cltq
	salq	$3, %rax
	testq	%rax, %rax
	je	.L869
	movq	PetscTrMalloc(%rip), %rax
	movq	-200(%rbp), %rdx
	addq	$240, %rdx
	movl	-44(%rbp), %ecx
	addl	$1, %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,8), %rbx
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1856, %esi
	movq	%rbx, %rdi
	call	*%rax
	jmp	.L870
.L869:
	movq	-200(%rbp), %rax
	movq	$0, 240(%rax)
	movl	$0, %eax
.L870:
	movl	%eax, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L868
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1856, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L868:
	.loc 4 1858 0
	movq	-200(%rbp), %rax
	movq	240(%rax), %rax
	movq	%rax, -40(%rbp)
	.loc 4 1859 0
	movq	-208(%rbp), %rax
	movq	%rax, -24(%rbp)
	.loc 4 1860 0
	movl	$0, -116(%rbp)
	jmp	.L871
.L879:
	.loc 4 1861 0
	movq	-104(%rbp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	-104(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -80(%rbp)
	addq	$4, -104(%rbp)
	.loc 4 1862 0
	movl	-80(%rbp), %eax
	imull	-88(%rbp), %eax
	movl	%eax, -48(%rbp)
	.loc 4 1863 0
	movl	-48(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L872
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1863, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L872:
	.loc 4 1864 0
	cmpl	$0, -52(%rbp)
	je	.L873
	.loc 4 1865 0
	movq	-208(%rbp), %rdx
	movl	-116(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-64(%rbp), %rax
	movl	(%rax), %eax
	imull	-88(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -24(%rbp)
.L873:
.LBB36:
	.loc 4 1867 0
	movabsq	$4607182418800017408, %rax
	movq	%rax, -264(%rbp)
	movl	$1, -268(%rbp)
	movl	-88(%rbp), %eax
	movl	%eax, -272(%rbp)
	movl	-48(%rbp), %eax
	movl	%eax, -276(%rbp)
	leaq	-272(%rbp), %rdi
	movq	-136(%rbp), %rsi
	leaq	-264(%rbp), %rcx
	leaq	-276(%rbp), %rdx
	leaq	-272(%rbp), %rax
	leaq	-268(%rbp), %rbx
	movq	%rbx, 32(%rsp)
	movq	-40(%rbp), %rbx
	movq	%rbx, 24(%rsp)
	leaq	-264(%rbp), %rbx
	movq	%rbx, 16(%rsp)
	leaq	-268(%rbp), %rbx
	movq	%rbx, 8(%rsp)
	movq	-24(%rbp), %rbx
	movq	%rbx, (%rsp)
	movq	%rdi, %r9
	movq	%rsi, %r8
	movq	%rax, %rsi
	movl	$.LC32, %edi
	call	dgemv_
.LBE36:
	.loc 4 1869 0
	movl	-80(%rbp), %eax
	imull	-76(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -136(%rbp)
	.loc 4 1870 0
	cmpl	$0, -52(%rbp)
	jne	.L874
	movl	-88(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -24(%rbp)
.L874:
	.loc 4 1871 0
	movq	-40(%rbp), %rax
	movq	%rax, -32(%rbp)
	.loc 4 1872 0
	movl	$0, -84(%rbp)
	jmp	.L875
.L878:
	.loc 4 1873 0
	movq	-216(%rbp), %rdx
	movq	-112(%rbp), %rax
	movl	(%rax), %eax
	imull	-88(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -192(%rbp)
	addq	$4, -112(%rbp)
	.loc 4 1874 0
	movl	$0, -44(%rbp)
	jmp	.L876
.L877:
	movl	-44(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-192(%rbp), %rax
	movl	-44(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	-192(%rbp), %rdx
	vmovsd	(%rdx), %xmm1
	movl	-44(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	-32(%rbp), %rdx
	vmovsd	(%rdx), %xmm0
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, (%rax)
	addl	$1, -44(%rbp)
.L876:
	movl	-44(%rbp), %eax
	cmpl	-88(%rbp), %eax
	jl	.L877
	.loc 4 1875 0
	movl	-88(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -32(%rbp)
	.loc 4 1872 0
	addl	$1, -84(%rbp)
.L875:
	movl	-84(%rbp), %eax
	cmpl	-80(%rbp), %eax
	jl	.L878
	.loc 4 1860 0
	addl	$1, -116(%rbp)
.L871:
	movl	-116(%rbp), %eax
	cmpl	-120(%rbp), %eax
	jl	.L879
.L843:
.LBE35:
	.loc 4 1880 0
	leaq	-208(%rbp), %rdx
	movq	-304(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L880
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1880, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L880:
	.loc 4 1881 0
	leaq	-216(%rbp), %rdx
	movq	-320(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L881
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1881, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L881:
	.loc 4 1882 0
	movq	-200(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vaddsd	%xmm0, %xmm0, %xmm1
	movq	-200(%rbp), %rax
	movl	224(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -124(%rbp)
	cmpl	$0, -124(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L882
	movl	-124(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.33932, %edx
	movl	$1882, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L825
.L882:
	.loc 4 1883 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L883
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L883
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L883:
	movl	$0, %eax
.L825:
	.loc 4 1884 0
	addq	$360, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE89:
	.size	MatMultTransposeAdd_SeqBAIJ, .-MatMultTransposeAdd_SeqBAIJ
.globl MatScale_SeqBAIJ
	.type	MatScale_SeqBAIJ, @function
MatScale_SeqBAIJ:
.LFB90:
	.loc 4 1889 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$72, %rsp
	movq	%rdi, -56(%rbp)
	vmovsd	%xmm0, -64(%rbp)
	.loc 4 1890 0
	movq	-56(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -32(%rbp)
	.loc 4 1891 0
	movq	-32(%rbp), %rax
	movl	224(%rax), %edx
	movq	-32(%rbp), %rax
	movl	128(%rax), %eax
	imull	%edx, %eax
	movl	%eax, -24(%rbp)
	.loc 4 1892 0
	movq	-64(%rbp), %rax
	movq	%rax, -40(%rbp)
	.loc 4 1894 0
	movl	$1, -44(%rbp)
	movl	-24(%rbp), %eax
	movl	%eax, -48(%rbp)
	.loc 4 1896 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L895
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L896
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.34647, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1896, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L894
.L895:
	nop
	jmp	.L894
.L896:
	nop
.L894:
	.loc 4 1897 0
	movq	-32(%rbp), %rax
	movq	168(%rax), %rdx
	leaq	-44(%rbp), %rcx
	leaq	-40(%rbp), %rbx
	leaq	-48(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	dscal_
	.loc 4 1898 0
	vcvtsi2sd	-24(%rbp), %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L890
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34647, %edx
	movl	$1898, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L891
.L890:
	.loc 4 1899 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L892
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L892
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L892:
	movl	$0, %eax
.L891:
	.loc 4 1900 0
	addq	$72, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE90:
	.size	MatScale_SeqBAIJ, .-MatScale_SeqBAIJ
	.section	.rodata
.LC33:
	.string	"No support for this norm yet"
	.text
.globl MatNorm_SeqBAIJ
	.type	MatNorm_SeqBAIJ, @function
MatNorm_SeqBAIJ:
.LFB91:
	.loc 4 1905 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$152, %rsp
	movq	%rdi, -104(%rbp)
	movl	%esi, -108(%rbp)
	movq	%rdx, -120(%rbp)
	.loc 4 1907 0
	movq	-104(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 1908 0
	movq	-80(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1909 0
	movl	$0, %eax
	movq	%rax, -64(%rbp)
	.loc 4 1910 0
	movq	-104(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -40(%rbp)
	movq	-80(%rbp), %rax
	movl	128(%rax), %eax
	movl	%eax, -36(%rbp)
	movq	-80(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -32(%rbp)
	.loc 4 1912 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L935
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L936
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.34724, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1912, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L934
.L935:
	nop
	jmp	.L934
.L936:
	nop
.L934:
	.loc 4 1913 0
	cmpl	$2, -108(%rbp)
	jne	.L899
	.loc 4 1914 0
	movl	$0, -52(%rbp)
	jmp	.L900
.L901:
	.loc 4 1918 0
	movq	-72(%rbp), %rax
	vmovsd	(%rax), %xmm1
	movq	-72(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	-64(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -64(%rbp)
	addq	$8, -72(%rbp)
	.loc 4 1914 0
	addl	$1, -52(%rbp)
.L900:
	movl	-32(%rbp), %eax
	imull	-36(%rbp), %eax
	cmpl	-52(%rbp), %eax
	jg	.L901
	.loc 4 1921 0
	vsqrtsd	-64(%rbp), %xmm0, %xmm0
	vucomisd	%xmm0, %xmm0
	jp	.L937
	je	.L902
.L937:
	vmovsd	-64(%rbp), %xmm0
	call	sqrt
.L902:
	vmovsd	%xmm0, -136(%rbp)
	movq	-136(%rbp), %rdx
	movq	-120(%rbp), %rax
	movq	%rdx, (%rax)
	jmp	.L903
.L899:
	.loc 4 1922 0
	cmpl	$0, -108(%rbp)
	jne	.L904
.LBB37:
	.loc 4 1924 0
	movq	-80(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -24(%rbp)
	.loc 4 1925 0
	movq	-104(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	addl	$1, %eax
	cltq
	salq	$3, %rax
	testq	%rax, %rax
	je	.L905
	movq	PetscTrMalloc(%rip), %rbx
	leaq	-96(%rbp), %rdx
	movq	-104(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	addl	$1, %eax
	cltq
	salq	$3, %rax
	movq	%rdx, %r9
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34724, %edx
	movl	$1925, %esi
	movq	%rax, %rdi
	call	*%rbx
	jmp	.L906
.L905:
	movq	$0, -96(%rbp)
	movl	$0, %eax
.L906:
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L907
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34724, %edx
	movl	$1925, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L908
.L907:
	.loc 4 1926 0
	movq	-104(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-96(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L909
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34724, %edx
	movl	$1926, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L908
.L909:
	.loc 4 1927 0
	movl	$0, -52(%rbp)
	jmp	.L910
.L915:
	.loc 4 1928 0
	movl	$0, -48(%rbp)
	jmp	.L911
.L914:
	.loc 4 1929 0
	movq	-24(%rbp), %rax
	movl	(%rax), %eax
	imull	-40(%rbp), %eax
	addl	-48(%rbp), %eax
	movl	%eax, -28(%rbp)
	.loc 4 1930 0
	movl	$0, -44(%rbp)
	jmp	.L912
.L913:
	.loc 4 1931 0
	movq	-96(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rbx
	movq	-96(%rbp), %rax
	movl	-28(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	vmovsd	(%rax), %xmm0
	vmovsd	%xmm0, -128(%rbp)
	movq	-72(%rbp), %rax
	vmovsd	(%rax), %xmm0
	call	PetscAbsScalar
	vaddsd	-128(%rbp), %xmm0, %xmm0
	vmovsd	%xmm0, (%rbx)
	addq	$8, -72(%rbp)
	.loc 4 1930 0
	addl	$1, -44(%rbp)
.L912:
	movl	-44(%rbp), %eax
	cmpl	-40(%rbp), %eax
	jl	.L913
	.loc 4 1928 0
	addl	$1, -48(%rbp)
.L911:
	movl	-48(%rbp), %eax
	cmpl	-40(%rbp), %eax
	jl	.L914
	.loc 4 1934 0
	addq	$4, -24(%rbp)
	.loc 4 1927 0
	addl	$1, -52(%rbp)
.L910:
	movl	-52(%rbp), %eax
	cmpl	-36(%rbp), %eax
	jl	.L915
	.loc 4 1936 0
	movq	-120(%rbp), %rax
	movl	$0, %edx
	movq	%rdx, (%rax)
	.loc 4 1937 0
	movl	$0, -48(%rbp)
	jmp	.L916
.L918:
	.loc 4 1938 0
	movq	-96(%rbp), %rax
	movl	-48(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	vmovsd	(%rax), %xmm0
	movq	-120(%rbp), %rax
	vmovsd	(%rax), %xmm1
	vucomisd	%xmm1, %xmm0
	seta	%al
	testb	%al, %al
	je	.L917
	movq	-96(%rbp), %rax
	movl	-48(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rdx
	movq	-120(%rbp), %rax
	movq	%rdx, (%rax)
.L917:
	.loc 4 1937 0
	addl	$1, -48(%rbp)
.L916:
	movq	-104(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	cmpl	-48(%rbp), %eax
	jg	.L918
	.loc 4 1940 0
	movq	-96(%rbp), %rax
	testq	%rax, %rax
	je	.L919
	movq	PetscTrFree(%rip), %rbx
	movq	-96(%rbp), %rax
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34724, %edx
	movl	$1940, %esi
	movq	%rax, %rdi
	call	*%rbx
	testl	%eax, %eax
	jne	.L920
	movq	$0, -96(%rbp)
	jmp	.L919
.L920:
	movl	$1, %eax
	jmp	.L921
.L919:
	movl	$0, %eax
.L921:
	movl	%eax, -84(%rbp)
	cmpl	$0, -84(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L903
	movl	-84(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34724, %edx
	movl	$1940, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L908
.L904:
.LBE37:
	.loc 4 1941 0
	cmpl	$3, -108(%rbp)
	jne	.L922
	.loc 4 1942 0
	movq	-120(%rbp), %rax
	movl	$0, %edx
	movq	%rdx, (%rax)
	.loc 4 1943 0
	movl	$0, -44(%rbp)
	jmp	.L923
.L931:
	.loc 4 1944 0
	movl	$0, -48(%rbp)
	jmp	.L924
.L930:
	.loc 4 1945 0
	movq	-80(%rbp), %rax
	movq	168(%rax), %rdx
	movq	-80(%rbp), %rax
	movq	136(%rax), %rax
	movl	-48(%rbp), %ecx
	movslq	%ecx, %rcx
	salq	$2, %rcx
	addq	%rcx, %rax
	movl	(%rax), %eax
	imull	-32(%rbp), %eax
	movslq	%eax, %rcx
	movl	-44(%rbp), %eax
	cltq
	leaq	(%rcx,%rax), %rax
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -72(%rbp)
	.loc 4 1946 0
	movl	$0, %eax
	movq	%rax, -64(%rbp)
	.loc 4 1947 0
	movl	$0, -52(%rbp)
	jmp	.L925
.L928:
	.loc 4 1948 0
	movl	$0, -28(%rbp)
	jmp	.L926
.L927:
	.loc 4 1949 0
	movq	-72(%rbp), %rax
	vmovsd	(%rax), %xmm0
	call	PetscAbsScalar
	vmovsd	-64(%rbp), %xmm1
	vaddsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -64(%rbp)
	.loc 4 1950 0
	movl	-40(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -72(%rbp)
	.loc 4 1948 0
	addl	$1, -28(%rbp)
.L926:
	movl	-28(%rbp), %eax
	cmpl	-40(%rbp), %eax
	jl	.L927
	.loc 4 1947 0
	addl	$1, -52(%rbp)
.L925:
	movq	-80(%rbp), %rax
	movq	136(%rax), %rax
	movl	-48(%rbp), %edx
	movslq	%edx, %rdx
	addq	$1, %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	movl	(%rax), %edx
	movq	-80(%rbp), %rax
	movq	136(%rax), %rax
	movl	-48(%rbp), %ecx
	movslq	%ecx, %rcx
	salq	$2, %rcx
	addq	%rcx, %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	cmpl	-52(%rbp), %eax
	jg	.L928
	.loc 4 1953 0
	movq	-120(%rbp), %rax
	vmovsd	(%rax), %xmm1
	vmovsd	-64(%rbp), %xmm0
	vucomisd	%xmm1, %xmm0
	seta	%al
	testb	%al, %al
	je	.L929
	movq	-120(%rbp), %rax
	movq	-64(%rbp), %rdx
	movq	%rdx, (%rax)
.L929:
	.loc 4 1944 0
	addl	$1, -48(%rbp)
.L924:
	movq	-80(%rbp), %rax
	movl	228(%rax), %eax
	cmpl	-48(%rbp), %eax
	jg	.L930
	.loc 4 1943 0
	addl	$1, -44(%rbp)
.L923:
	movl	-44(%rbp), %eax
	cmpl	-40(%rbp), %eax
	jl	.L931
	jmp	.L903
.L922:
	.loc 4 1956 0
	movq	$.LC33, 8(%rsp)
	movl	$0, (%rsp)
	movl	$56, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34724, %edx
	movl	$1956, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L908
.L903:
	.loc 4 1957 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L932
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L932
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L932:
	movl	$0, %eax
.L908:
	.loc 4 1958 0
	addq	$152, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE91:
	.size	MatNorm_SeqBAIJ, .-MatNorm_SeqBAIJ
.globl MatEqual_SeqBAIJ
	.type	MatEqual_SeqBAIJ, @function
MatEqual_SeqBAIJ:
.LFB92:
	.loc 4 1964 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$88, %rsp
	movq	%rdi, -56(%rbp)
	movq	%rsi, -64(%rbp)
	movq	%rdx, -72(%rbp)
	.loc 4 1965 0
	movq	-56(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -40(%rbp)
	movq	-64(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -32(%rbp)
	.loc 4 1968 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L954
	.cfi_offset 3, -24
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L955
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.34937, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$1968, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L953
.L954:
	nop
	jmp	.L953
.L955:
	nop
.L953:
	.loc 4 1970 0
	movq	-56(%rbp), %rax
	movq	456(%rax), %rax
	movl	8(%rax), %edx
	movq	-64(%rbp), %rax
	movq	456(%rax), %rax
	movl	8(%rax), %eax
	cmpl	%eax, %edx
	jne	.L940
	movq	-56(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %edx
	movq	-64(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	cmpl	%eax, %edx
	jne	.L940
	movq	-56(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %edx
	movq	-64(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	cmpl	%eax, %edx
	jne	.L940
	movq	-40(%rbp), %rax
	movl	128(%rax), %edx
	movq	-32(%rbp), %rax
	movl	128(%rax), %eax
	cmpl	%eax, %edx
	je	.L941
.L940:
	.loc 4 1971 0
	movq	-72(%rbp), %rax
	movl	$0, (%rax)
	.loc 4 1972 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L942
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L942
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L942:
	movl	$0, %eax
	jmp	.L943
.L941:
	.loc 4 1976 0
	movq	-40(%rbp), %rax
	movl	228(%rax), %eax
	addl	$1, %eax
	cltq
	leaq	0(,%rax,4), %rsi
	movq	-32(%rbp), %rax
	movq	136(%rax), %rbx
	movq	-40(%rbp), %rax
	movq	136(%rax), %rax
	movq	-72(%rbp), %rdx
	movq	%rdx, %rcx
	movq	%rsi, %rdx
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	PetscMemcmp
	movl	%eax, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L944
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34937, %edx
	movl	$1976, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L943
.L944:
	.loc 4 1977 0
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	testl	%eax, %eax
	jne	.L945
	.loc 4 1978 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L946
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L946
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L946:
	movl	$0, %eax
	jmp	.L943
.L945:
	.loc 4 1982 0
	movq	-40(%rbp), %rax
	movl	128(%rax), %eax
	cltq
	leaq	0(,%rax,4), %rsi
	movq	-32(%rbp), %rax
	movq	144(%rax), %rbx
	movq	-40(%rbp), %rax
	movq	144(%rax), %rax
	movq	-72(%rbp), %rdx
	movq	%rdx, %rcx
	movq	%rsi, %rdx
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	PetscMemcmp
	movl	%eax, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L947
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34937, %edx
	movl	$1982, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L943
.L947:
	.loc 4 1983 0
	movq	-72(%rbp), %rax
	movl	(%rax), %eax
	testl	%eax, %eax
	jne	.L948
	.loc 4 1984 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L949
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L949
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L949:
	movl	$0, %eax
	jmp	.L943
.L948:
	.loc 4 1987 0
	movq	-40(%rbp), %rax
	movl	128(%rax), %edx
	movq	-56(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	imull	%eax, %edx
	movq	-64(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	imull	%edx, %eax
	cltq
	leaq	0(,%rax,8), %rsi
	movq	-32(%rbp), %rax
	movq	168(%rax), %rbx
	movq	-40(%rbp), %rax
	movq	168(%rax), %rax
	movq	-72(%rbp), %rdx
	movq	%rdx, %rcx
	movq	%rsi, %rdx
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	PetscMemcmp
	movl	%eax, -20(%rbp)
	cmpl	$0, -20(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L950
	movl	-20(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.34937, %edx
	movl	$1987, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L943
.L950:
	.loc 4 1988 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L951
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L951
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L951:
	movl	$0, %eax
.L943:
	.loc 4 1990 0
	addq	$88, %rsp
	popq	%rbx
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE92:
	.size	MatEqual_SeqBAIJ, .-MatEqual_SeqBAIJ
	.section	.rodata
.LC34:
	.string	"Not for factored matrix"
	.align 8
.LC35:
	.string	"Nonconforming matrix and vector"
	.text
.globl MatGetDiagonal_SeqBAIJ
	.type	MatGetDiagonal_SeqBAIJ, @function
MatGetDiagonal_SeqBAIJ:
.LFB93:
	.loc 4 1995 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	addq	$-128, %rsp
	movq	%rdi, -104(%rbp)
	movq	%rsi, -112(%rbp)
	.loc 4 1996 0
	movq	-104(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 1999 0
	movl	$0, %eax
	movq	%rax, -24(%rbp)
	.loc 4 2002 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L976
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L977
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.35139, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$2002, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L975
.L976:
	nop
	jmp	.L975
.L977:
	nop
.L975:
	.loc 4 2003 0
	movq	-104(%rbp), %rax
	movl	480(%rax), %eax
	testl	%eax, %eax
	je	.L958
	movq	$.LC34, 8(%rsp)
	movl	$0, (%rsp)
	movl	$73, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35139, %edx
	movl	$2003, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L959
.L958:
	.loc 4 2004 0
	movq	-104(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -52(%rbp)
	.loc 4 2005 0
	movq	-80(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -16(%rbp)
	.loc 4 2006 0
	movq	-80(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -48(%rbp)
	.loc 4 2007 0
	movq	-80(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -40(%rbp)
	.loc 4 2008 0
	movq	-80(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -32(%rbp)
	.loc 4 2009 0
	movq	-80(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -28(%rbp)
	.loc 4 2011 0
	vmovsd	-24(%rbp), %xmm0
	movq	-112(%rbp), %rax
	movq	%rax, %rdi
	call	VecSet
	movl	%eax, -72(%rbp)
	cmpl	$0, -72(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L960
	movl	-72(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35139, %edx
	movl	$2011, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L959
.L960:
	.loc 4 2012 0
	leaq	-96(%rbp), %rdx
	movq	-112(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArray
	movl	%eax, -72(%rbp)
	cmpl	$0, -72(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L961
	movl	-72(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35139, %edx
	movl	$2012, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L959
.L961:
	.loc 4 2013 0
	leaq	-84(%rbp), %rdx
	movq	-112(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetLocalSize
	movl	%eax, -72(%rbp)
	cmpl	$0, -72(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L962
	movl	-72(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35139, %edx
	movl	$2013, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L959
.L962:
	.loc 4 2014 0
	movq	-104(%rbp), %rax
	movq	456(%rax), %rax
	movl	8(%rax), %edx
	movl	-84(%rbp), %eax
	cmpl	%eax, %edx
	je	.L963
	movq	$.LC35, 8(%rsp)
	movl	$0, (%rsp)
	movl	$60, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35139, %edx
	movl	$2014, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L959
.L963:
	.loc 4 2015 0
	movl	$0, -68(%rbp)
	jmp	.L964
.L971:
	.loc 4 2016 0
	movl	-68(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -64(%rbp)
	jmp	.L965
.L970:
	.loc 4 2017 0
	movl	-64(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-40(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-68(%rbp), %eax
	jne	.L966
	.loc 4 2018 0
	movl	-68(%rbp), %eax
	imull	-52(%rbp), %eax
	movl	%eax, -56(%rbp)
	.loc 4 2019 0
	movl	-64(%rbp), %eax
	imull	-28(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-16(%rbp), %rax
	movq	%rax, -8(%rbp)
	.loc 4 2020 0
	movl	$0, -60(%rbp)
	jmp	.L967
.L968:
	movq	-96(%rbp), %rax
	movl	-56(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rdx
	movl	-60(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, (%rdx)
	movl	-52(%rbp), %eax
	addl	$1, %eax
	addl	%eax, -60(%rbp)
	addl	$1, -56(%rbp)
.L967:
	movl	-60(%rbp), %eax
	cmpl	-28(%rbp), %eax
	jl	.L968
	.loc 4 2021 0
	jmp	.L969
.L966:
	.loc 4 2016 0
	addl	$1, -64(%rbp)
.L965:
	movl	-68(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-48(%rbp), %rax
	movl	(%rax), %eax
	cmpl	-64(%rbp), %eax
	jg	.L970
.L969:
	.loc 4 2015 0
	addl	$1, -68(%rbp)
.L964:
	movl	-68(%rbp), %eax
	cmpl	-32(%rbp), %eax
	jl	.L971
	.loc 4 2025 0
	leaq	-96(%rbp), %rdx
	movq	-112(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArray
	movl	%eax, -72(%rbp)
	cmpl	$0, -72(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L972
	movl	-72(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35139, %edx
	movl	$2025, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L959
.L972:
	.loc 4 2026 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L973
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L973
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L973:
	movl	$0, %eax
.L959:
	.loc 4 2027 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE93:
	.size	MatGetDiagonal_SeqBAIJ, .-MatGetDiagonal_SeqBAIJ
	.section	.rodata
	.align 8
.LC36:
	.string	"Left scaling vector wrong length"
	.align 8
.LC37:
	.string	"Right scaling vector wrong length"
	.text
.globl MatDiagonalScale_SeqBAIJ
	.type	MatDiagonalScale_SeqBAIJ, @function
MatDiagonalScale_SeqBAIJ:
.LFB94:
	.loc 4 2032 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$192, %rsp
	movq	%rdi, -152(%rbp)
	movq	%rsi, -160(%rbp)
	movq	%rdx, -168(%rbp)
	.loc 4 2033 0
	movq	-152(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -112(%rbp)
	.loc 4 2041 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L1008
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L1009
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.35288, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$2041, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L1007
.L1008:
	nop
	jmp	.L1007
.L1009:
	nop
.L1007:
	.loc 4 2042 0
	movq	-112(%rbp), %rax
	movq	136(%rax), %rax
	movq	%rax, -16(%rbp)
	.loc 4 2043 0
	movq	-112(%rbp), %rax
	movq	144(%rax), %rax
	movq	%rax, -8(%rbp)
	.loc 4 2044 0
	movq	-112(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rax, -80(%rbp)
	.loc 4 2045 0
	movq	-152(%rbp), %rax
	movq	456(%rax), %rax
	movl	4(%rax), %eax
	movl	%eax, -44(%rbp)
	.loc 4 2046 0
	movq	-152(%rbp), %rax
	movq	464(%rax), %rax
	movl	4(%rax), %eax
	movl	%eax, -40(%rbp)
	.loc 4 2047 0
	movq	-152(%rbp), %rax
	movq	456(%rax), %rax
	movl	32(%rax), %eax
	movl	%eax, -28(%rbp)
	.loc 4 2048 0
	movq	-112(%rbp), %rax
	movl	228(%rax), %eax
	movl	%eax, -36(%rbp)
	.loc 4 2049 0
	movq	-112(%rbp), %rax
	movl	224(%rax), %eax
	movl	%eax, -24(%rbp)
	.loc 4 2050 0
	cmpq	$0, -160(%rbp)
	je	.L980
	.loc 4 2051 0
	leaq	-120(%rbp), %rdx
	movq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L981
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2051, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L981:
	.loc 4 2052 0
	leaq	-132(%rbp), %rdx
	movq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetLocalSize
	movl	%eax, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L983
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2052, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L983:
	.loc 4 2053 0
	movl	-132(%rbp), %eax
	cmpl	-44(%rbp), %eax
	je	.L984
	movq	$.LC36, 8(%rsp)
	movl	$0, (%rsp)
	movl	$60, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2053, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L984:
	.loc 4 2054 0
	movl	$0, -60(%rbp)
	jmp	.L985
.L990:
	.loc 4 2055 0
	movl	-60(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-16(%rbp), %rax
	movl	(%rax), %edx
	movl	-60(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-16(%rbp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	movl	%eax, -48(%rbp)
	.loc 4 2056 0
	movq	-120(%rbp), %rdx
	movl	-60(%rbp), %eax
	imull	-28(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -104(%rbp)
	.loc 4 2057 0
	movl	-60(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-16(%rbp), %rax
	movl	(%rax), %eax
	imull	-24(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-80(%rbp), %rax
	movq	%rax, -72(%rbp)
	.loc 4 2058 0
	movl	$0, -56(%rbp)
	jmp	.L986
.L989:
	.loc 4 2059 0
	movl	$0, -52(%rbp)
	jmp	.L987
.L988:
	.loc 4 2060 0
	movq	-72(%rbp), %rcx
	vmovsd	(%rcx), %xmm1
	movl	-52(%rbp), %eax
	movl	%eax, %edx
	sarl	$31, %edx
	idivl	-28(%rbp)
	movl	%edx, %eax
	cltq
	salq	$3, %rax
	addq	-104(%rbp), %rax
	vmovsd	(%rax), %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, (%rcx)
	addq	$8, -72(%rbp)
	.loc 4 2059 0
	addl	$1, -52(%rbp)
.L987:
	movl	-52(%rbp), %eax
	cmpl	-24(%rbp), %eax
	jl	.L988
	.loc 4 2058 0
	addl	$1, -56(%rbp)
.L986:
	movl	-56(%rbp), %eax
	cmpl	-48(%rbp), %eax
	jl	.L989
	.loc 4 2054 0
	addl	$1, -60(%rbp)
.L985:
	movl	-60(%rbp), %eax
	cmpl	-36(%rbp), %eax
	jl	.L990
	.loc 4 2064 0
	leaq	-120(%rbp), %rdx
	movq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L991
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2064, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L991:
	.loc 4 2065 0
	movq	-112(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L980
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2065, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L980:
	.loc 4 2068 0
	cmpq	$0, -168(%rbp)
	je	.L992
	.loc 4 2069 0
	leaq	-128(%rbp), %rdx
	movq	-168(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetArrayRead
	movl	%eax, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L993
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2069, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L993:
	.loc 4 2070 0
	leaq	-136(%rbp), %rdx
	movq	-168(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecGetLocalSize
	movl	%eax, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L994
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2070, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L994:
	.loc 4 2071 0
	movl	-136(%rbp), %eax
	cmpl	-40(%rbp), %eax
	je	.L995
	movq	$.LC37, 8(%rsp)
	movl	$0, (%rsp)
	movl	$60, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2071, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L995:
	.loc 4 2072 0
	movl	$0, -60(%rbp)
	jmp	.L996
.L1003:
	.loc 4 2073 0
	movl	-60(%rbp), %eax
	cltq
	salq	$2, %rax
	addq	-16(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -20(%rbp)
	.loc 4 2074 0
	movl	-60(%rbp), %eax
	cltq
	addq	$1, %rax
	salq	$2, %rax
	addq	-16(%rbp), %rax
	movl	(%rax), %eax
	subl	-20(%rbp), %eax
	movl	%eax, -48(%rbp)
	.loc 4 2075 0
	movl	-24(%rbp), %eax
	imull	-20(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-80(%rbp), %rax
	movq	%rax, -72(%rbp)
	.loc 4 2076 0
	movl	$0, -56(%rbp)
	jmp	.L997
.L1002:
	.loc 4 2077 0
	movq	-128(%rbp), %rdx
	movl	-56(%rbp), %eax
	movl	-20(%rbp), %ecx
	leal	(%rcx,%rax), %eax
	cltq
	salq	$2, %rax
	addq	-8(%rbp), %rax
	movl	(%rax), %eax
	imull	-28(%rbp), %eax
	cltq
	salq	$3, %rax
	leaq	(%rdx,%rax), %rax
	movq	%rax, -96(%rbp)
	.loc 4 2078 0
	movl	$0, -52(%rbp)
	jmp	.L998
.L1001:
	.loc 4 2079 0
	movl	-52(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-96(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -88(%rbp)
	.loc 4 2080 0
	movl	$0, -32(%rbp)
	jmp	.L999
.L1000:
	movl	-32(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	-72(%rbp), %rax
	movl	-32(%rbp), %edx
	movslq	%edx, %rdx
	salq	$3, %rdx
	addq	-72(%rbp), %rdx
	vmovsd	(%rdx), %xmm0
	vmulsd	-88(%rbp), %xmm0, %xmm0
	vmovsd	%xmm0, (%rax)
	addl	$1, -32(%rbp)
.L999:
	movl	-32(%rbp), %eax
	cmpl	-28(%rbp), %eax
	jl	.L1000
	.loc 4 2081 0
	movl	-28(%rbp), %eax
	cltq
	salq	$3, %rax
	addq	%rax, -72(%rbp)
	.loc 4 2078 0
	addl	$1, -52(%rbp)
.L998:
	movl	-52(%rbp), %eax
	cmpl	-28(%rbp), %eax
	jl	.L1001
	.loc 4 2076 0
	addl	$1, -56(%rbp)
.L997:
	movl	-56(%rbp), %eax
	cmpl	-48(%rbp), %eax
	jl	.L1002
	.loc 4 2072 0
	addl	$1, -60(%rbp)
.L996:
	movl	-60(%rbp), %eax
	cmpl	-36(%rbp), %eax
	jl	.L1003
	.loc 4 2085 0
	leaq	-128(%rbp), %rdx
	movq	-168(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	VecRestoreArrayRead
	movl	%eax, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L1004
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2085, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L1004:
	.loc 4 2086 0
	movq	-112(%rbp), %rax
	movl	128(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	vmovsd	_TotalFlops(%rip), %xmm1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm0, _TotalFlops(%rip)
	movl	$0, -64(%rbp)
	cmpl	$0, -64(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L992
	movl	-64(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35288, %edx
	movl	$2086, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L982
.L992:
	.loc 4 2088 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L1005
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L1005
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L1005:
	movl	$0, %eax
.L982:
	.loc 4 2089 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE94:
	.size	MatDiagonalScale_SeqBAIJ, .-MatDiagonalScale_SeqBAIJ
.globl MatGetInfo_SeqBAIJ
	.type	MatGetInfo_SeqBAIJ, @function
MatGetInfo_SeqBAIJ:
.LFB95:
	.loc 4 2095 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movl	%esi, -28(%rbp)
	movq	%rdx, -40(%rbp)
	.loc 4 2096 0
	movq	-24(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -8(%rbp)
	.loc 4 2098 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L1017
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L1018
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.35495, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$2098, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L1016
.L1017:
	nop
	jmp	.L1016
.L1018:
	nop
.L1016:
	.loc 4 2099 0
	movq	-8(%rbp), %rax
	movl	224(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	movq	-40(%rbp), %rax
	vmovsd	%xmm0, (%rax)
	.loc 4 2100 0
	movq	-8(%rbp), %rax
	movl	224(%rax), %edx
	movq	-8(%rbp), %rax
	movl	16(%rax), %eax
	imull	%edx, %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	movq	-40(%rbp), %rax
	vmovsd	%xmm0, 8(%rax)
	.loc 4 2101 0
	movq	-8(%rbp), %rax
	movl	224(%rax), %edx
	movq	-8(%rbp), %rax
	movl	128(%rax), %eax
	imull	%edx, %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	movq	-40(%rbp), %rax
	vmovsd	%xmm0, 16(%rax)
	.loc 4 2102 0
	movq	-40(%rbp), %rax
	vmovsd	8(%rax), %xmm0
	movq	-40(%rbp), %rax
	vmovsd	16(%rax), %xmm1
	vsubsd	%xmm1, %xmm0, %xmm0
	movq	-40(%rbp), %rax
	vmovsd	%xmm0, 24(%rax)
	.loc 4 2103 0
	movq	-24(%rbp), %rax
	movl	492(%rax), %eax
	vcvtsi2sd	%eax, %xmm0, %xmm0
	movq	-40(%rbp), %rax
	vmovsd	%xmm0, 40(%rax)
	.loc 4 2104 0
	movq	-24(%rbp), %rax
	movq	552(%rax), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, 48(%rax)
	.loc 4 2105 0
	movq	-24(%rbp), %rax
	movq	40(%rax), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, 32(%rax)
	.loc 4 2106 0
	movq	-24(%rbp), %rax
	movl	480(%rax), %eax
	testl	%eax, %eax
	je	.L1012
	.loc 4 2107 0
	movq	-24(%rbp), %rax
	movq	560(%rax), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, 56(%rax)
	.loc 4 2108 0
	movq	-24(%rbp), %rax
	movq	568(%rax), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, 64(%rax)
	.loc 4 2109 0
	movq	-24(%rbp), %rax
	movq	576(%rax), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, 72(%rax)
	jmp	.L1013
.L1012:
	.loc 4 2111 0
	movq	-40(%rbp), %rax
	movl	$0, %edx
	movq	%rdx, 56(%rax)
	.loc 4 2112 0
	movq	-40(%rbp), %rax
	movl	$0, %edx
	movq	%rdx, 64(%rax)
	.loc 4 2113 0
	movq	-40(%rbp), %rax
	movl	$0, %edx
	movq	%rdx, 72(%rax)
.L1013:
	.loc 4 2115 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L1014
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L1014
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L1014:
	movl	$0, %eax
	.loc 4 2116 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE95:
	.size	MatGetInfo_SeqBAIJ, .-MatGetInfo_SeqBAIJ
.globl MatZeroEntries_SeqBAIJ
	.type	MatZeroEntries_SeqBAIJ, @function
MatZeroEntries_SeqBAIJ:
.LFB96:
	.loc 4 2122 0
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	.loc 4 2123 0
	movq	-24(%rbp), %rax
	movq	472(%rax), %rax
	movq	%rax, -16(%rbp)
	.loc 4 2126 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L1026
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	cmpl	$63, %eax
	jg	.L1027
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$__func__.35575, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$.LC6, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$.LC3, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$2126, (%rax,%rdx,4)
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	addl	$1, %edx
	movl	%edx, 1792(%rax)
	jmp	.L1025
.L1026:
	nop
	jmp	.L1025
.L1027:
	nop
.L1025:
	.loc 4 2127 0
	movq	-16(%rbp), %rax
	movl	224(%rax), %edx
	movq	-16(%rbp), %rax
	movq	136(%rax), %rcx
	movq	-16(%rbp), %rax
	movl	228(%rax), %eax
	cltq
	salq	$2, %rax
	leaq	(%rcx,%rax), %rax
	movl	(%rax), %eax
	imull	%edx, %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-16(%rbp), %rax
	movq	168(%rax), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	PetscMemzero
	movl	%eax, -4(%rbp)
	cmpl	$0, -4(%rbp)
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L1021
	movl	-4(%rbp), %eax
	movq	$.LC5, 8(%rsp)
	movl	$1, (%rsp)
	movl	%eax, %r9d
	movl	$.LC3, %r8d
	movl	$.LC6, %ecx
	movl	$__func__.35575, %edx
	movl	$2127, %esi
	movl	$1140850689, %edi
	movl	$0, %eax
	call	PetscError
	jmp	.L1022
.L1021:
	.loc 4 2128 0
	movq	petscstack(%rip), %rax
	testq	%rax, %rax
	je	.L1023
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %eax
	testl	%eax, %eax
	jle	.L1023
	movq	petscstack(%rip), %rax
	movl	1792(%rax), %edx
	subl	$1, %edx
	movl	%edx, 1792(%rax)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$64, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	subq	$-128, %rdx
	movq	$0, (%rax,%rdx,8)
	movq	petscstack(%rip), %rax
	movq	petscstack(%rip), %rdx
	movl	1792(%rdx), %edx
	movslq	%edx, %rdx
	addq	$384, %rdx
	movl	$0, (%rax,%rdx,4)
.L1023:
	movl	$0, %eax
.L1022:
	.loc 4 2129 0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE96:
	.size	MatZeroEntries_SeqBAIJ, .-MatZeroEntries_SeqBAIJ
	.section	.rodata
	.align 16
	.type	__func__.35575, @object
	.size	__func__.35575, 23
__func__.35575:
	.string	"MatZeroEntries_SeqBAIJ"
	.align 16
	.type	__func__.35495, @object
	.size	__func__.35495, 19
__func__.35495:
	.string	"MatGetInfo_SeqBAIJ"
	.align 16
	.type	__func__.35288, @object
	.size	__func__.35288, 25
__func__.35288:
	.string	"MatDiagonalScale_SeqBAIJ"
	.align 16
	.type	__func__.15175, @object
	.size	__func__.15175, 16
__func__.15175:
	.string	"VecGetArrayRead"
	.align 16
	.type	__func__.15257, @object
	.size	__func__.15257, 20
__func__.15257:
	.string	"VecRestoreArrayRead"
	.align 16
	.type	__func__.35139, @object
	.size	__func__.35139, 23
__func__.35139:
	.string	"MatGetDiagonal_SeqBAIJ"
	.type	__func__.15326, @object
	.size	__func__.15326, 12
__func__.15326:
	.string	"VecGetArray"
	.align 16
	.type	__func__.15408, @object
	.size	__func__.15408, 16
__func__.15408:
	.string	"VecRestoreArray"
	.align 16
	.type	__func__.34937, @object
	.size	__func__.34937, 17
__func__.34937:
	.string	"MatEqual_SeqBAIJ"
	.align 16
	.type	__func__.34724, @object
	.size	__func__.34724, 16
__func__.34724:
	.string	"MatNorm_SeqBAIJ"
	.align 16
	.type	__func__.34647, @object
	.size	__func__.34647, 17
__func__.34647:
	.string	"MatScale_SeqBAIJ"
	.align 16
	.type	__func__.33932, @object
	.size	__func__.33932, 28
__func__.33932:
	.string	"MatMultTransposeAdd_SeqBAIJ"
	.align 32
	.type	__func__.33195, @object
	.size	__func__.33195, 37
__func__.33195:
	.string	"MatMultHermitianTransposeAdd_SeqBAIJ"
	.align 16
	.type	__func__.33103, @object
	.size	__func__.33103, 25
__func__.33103:
	.string	"MatMultTranspose_SeqBAIJ"
	.align 32
	.type	__func__.33035, @object
	.size	__func__.33035, 34
__func__.33035:
	.string	"MatMultHermitianTranspose_SeqBAIJ"
	.align 16
	.type	__func__.32859, @object
	.size	__func__.32859, 21
__func__.32859:
	.string	"MatMultAdd_SeqBAIJ_N"
	.align 16
	.type	__func__.32456, @object
	.size	__func__.32456, 21
__func__.32456:
	.string	"MatMultAdd_SeqBAIJ_7"
	.type	__func__.12688, @object
	.size	__func__.12688, 12
__func__.12688:
	.string	"PetscMemcpy"
	.align 16
	.type	__func__.32097, @object
	.size	__func__.32097, 21
__func__.32097:
	.string	"MatMultAdd_SeqBAIJ_6"
	.align 16
	.type	__func__.31786, @object
	.size	__func__.31786, 21
__func__.31786:
	.string	"MatMultAdd_SeqBAIJ_5"
	.align 16
	.type	__func__.31515, @object
	.size	__func__.31515, 21
__func__.31515:
	.string	"MatMultAdd_SeqBAIJ_4"
	.align 16
	.type	__func__.31276, @object
	.size	__func__.31276, 21
__func__.31276:
	.string	"MatMultAdd_SeqBAIJ_3"
	.align 16
	.type	__func__.31048, @object
	.size	__func__.31048, 21
__func__.31048:
	.string	"MatMultAdd_SeqBAIJ_2"
	.align 16
	.type	__func__.30821, @object
	.size	__func__.30821, 21
__func__.30821:
	.string	"MatMultAdd_SeqBAIJ_1"
	.align 16
	.type	__func__.30632, @object
	.size	__func__.30632, 18
__func__.30632:
	.string	"MatMult_SeqBAIJ_N"
	.align 16
	.type	__func__.29563, @object
	.size	__func__.29563, 24
__func__.29563:
	.string	"MatMult_SeqBAIJ_15_ver4"
	.align 16
	.type	__func__.28485, @object
	.size	__func__.28485, 24
__func__.28485:
	.string	"MatMult_SeqBAIJ_15_ver3"
	.align 16
	.type	__func__.27446, @object
	.size	__func__.27446, 24
__func__.27446:
	.string	"MatMult_SeqBAIJ_15_ver2"
	.align 16
	.type	__func__.27211, @object
	.size	__func__.27211, 24
__func__.27211:
	.string	"MatMult_SeqBAIJ_15_ver1"
	.align 16
	.type	__func__.26835, @object
	.size	__func__.26835, 18
__func__.26835:
	.string	"MatMult_SeqBAIJ_7"
	.align 16
	.type	__func__.26516, @object
	.size	__func__.26516, 18
__func__.26516:
	.string	"MatMult_SeqBAIJ_6"
	.align 16
	.type	__func__.26244, @object
	.size	__func__.26244, 18
__func__.26244:
	.string	"MatMult_SeqBAIJ_5"
	.align 16
	.type	__func__.26011, @object
	.size	__func__.26011, 18
__func__.26011:
	.string	"MatMult_SeqBAIJ_4"
	.align 16
	.type	__func__.25809, @object
	.size	__func__.25809, 18
__func__.25809:
	.string	"MatMult_SeqBAIJ_3"
	.align 16
	.type	__func__.25630, @object
	.size	__func__.25630, 18
__func__.25630:
	.string	"MatMult_SeqBAIJ_2"
	.align 16
	.type	__func__.25437, @object
	.size	__func__.25437, 18
__func__.25437:
	.string	"MatMult_SeqBAIJ_1"
	.align 16
	.type	__func__.25328, @object
	.size	__func__.25328, 26
__func__.25328:
	.string	"MatGetSubMatrices_SeqBAIJ"
	.align 16
	.type	__func__.25037, @object
	.size	__func__.25037, 24
__func__.25037:
	.string	"MatGetSubMatrix_SeqBAIJ"
	.align 32
	.type	__func__.24627, @object
	.size	__func__.24627, 32
__func__.24627:
	.string	"MatGetSubMatrix_SeqBAIJ_Private"
	.align 16
	.type	__func__.24258, @object
	.size	__func__.24258, 27
__func__.24258:
	.string	"MatIncreaseOverlap_SeqBAIJ"
	.align 16
.LC1:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.align 8
.LC14:
	.long	0
	.long	1075838976
	.align 8
.LC15:
	.long	0
	.long	-1073741824
	.align 8
.LC16:
	.long	0
	.long	1077018624
	.align 8
.LC17:
	.long	0
	.long	-1073217536
	.align 8
.LC18:
	.long	0
	.long	1077936128
	.align 8
.LC19:
	.long	0
	.long	-1072693248
	.align 8
.LC20:
	.long	0
	.long	1078525952
	.align 8
.LC21:
	.long	0
	.long	-1072431104
	.align 8
.LC22:
	.long	0
	.long	1079115776
	.align 8
.LC23:
	.long	0
	.long	-1072168960
	.align 8
.LC24:
	.long	0
	.long	1079541760
	.align 8
.LC25:
	.long	0
	.long	-1071906816
	.align 8
.LC26:
	.long	0
	.long	1081876480
	.align 8
.LC27:
	.long	0
	.long	-1070727168
	.align 8
.LC30:
	.long	0
	.long	1074790400
	.text
.Letext0:
	.file 5 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
	.file 6 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h"
	.file 7 "/usr/include/bits/types.h"
	.file 8 "/usr/include/libio.h"
	.file 9 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h"
	.file 10 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"
	.file 11 "/usr/include/sys/types.h"
	.file 12 "/usr/include/stdint.h"
	.file 13 "/home/dpnkarthik/petsc-rnet/include/petscerror.h"
	.file 14 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
	.file 15 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
	.file 16 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
	.file 17 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
	.file 18 "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"
	.file 19 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/baij/seq/baij.h"
	.file 20 "/home/dpnkarthik/petsc-rnet/include/petscbt.h"
	.file 21 "/usr/include/stdio.h"
	.file 22 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"
	.file 23 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Fpublic.h"
	.section	.debug_info
	.long	0x7a9f
	.value	0x3
	.long	.Ldebug_abbrev0
	.byte	0x8
	.uleb128 0x1
	.long	.LASF748
	.byte	0x1
	.long	.LASF749
	.long	.LASF750
	.quad	.Ltext0
	.quad	.Letext0
	.long	.Ldebug_line0
	.uleb128 0x2
	.byte	0x4
	.byte	0x5
	.string	"int"
	.uleb128 0x3
	.long	.LASF0
	.byte	0x5
	.byte	0x7c
	.long	0x2d
	.uleb128 0x4
	.byte	0x8
	.long	0x34
	.uleb128 0x4
	.byte	0x8
	.long	0x2d
	.uleb128 0x5
	.byte	0x8
	.uleb128 0x6
	.long	.LASF1
	.byte	0x5
	.value	0x114
	.long	0x2d
	.uleb128 0x7
	.byte	0x8
	.byte	0x5
	.long	.LASF2
	.uleb128 0x7
	.byte	0x8
	.byte	0x5
	.long	.LASF3
	.uleb128 0x8
	.long	.LASF9
	.byte	0x14
	.byte	0x5
	.value	0x184
	.long	0xb6
	.uleb128 0x9
	.long	.LASF4
	.byte	0x5
	.value	0x185
	.long	0x2d
	.sleb128 0
	.uleb128 0x9
	.long	.LASF5
	.byte	0x5
	.value	0x186
	.long	0x2d
	.sleb128 4
	.uleb128 0x9
	.long	.LASF6
	.byte	0x5
	.value	0x187
	.long	0x2d
	.sleb128 8
	.uleb128 0x9
	.long	.LASF7
	.byte	0x5
	.value	0x188
	.long	0x2d
	.sleb128 12
	.uleb128 0x9
	.long	.LASF8
	.byte	0x5
	.value	0x189
	.long	0x2d
	.sleb128 16
	.byte	0x0
	.uleb128 0x6
	.long	.LASF9
	.byte	0x5
	.value	0x18b
	.long	0x67
	.uleb128 0x4
	.byte	0x8
	.long	0xb6
	.uleb128 0x4
	.byte	0x8
	.long	0x4b
	.uleb128 0x7
	.byte	0x8
	.byte	0x4
	.long	.LASF10
	.uleb128 0x3
	.long	.LASF11
	.byte	0x6
	.byte	0xd3
	.long	0xe0
	.uleb128 0x7
	.byte	0x8
	.byte	0x7
	.long	.LASF12
	.uleb128 0x7
	.byte	0x1
	.byte	0x8
	.long	.LASF13
	.uleb128 0x7
	.byte	0x2
	.byte	0x7
	.long	.LASF14
	.uleb128 0x7
	.byte	0x4
	.byte	0x7
	.long	.LASF15
	.uleb128 0x7
	.byte	0x1
	.byte	0x6
	.long	.LASF16
	.uleb128 0x7
	.byte	0x2
	.byte	0x5
	.long	.LASF17
	.uleb128 0x3
	.long	.LASF18
	.byte	0x7
	.byte	0x8d
	.long	0x59
	.uleb128 0x3
	.long	.LASF19
	.byte	0x7
	.byte	0x8e
	.long	0x59
	.uleb128 0x4
	.byte	0x8
	.long	0x126
	.uleb128 0x7
	.byte	0x1
	.byte	0x6
	.long	.LASF20
	.uleb128 0x8
	.long	.LASF21
	.byte	0xd8
	.byte	0x8
	.value	0x10f
	.long	0x2c9
	.uleb128 0x9
	.long	.LASF22
	.byte	0x8
	.value	0x110
	.long	0x2d
	.sleb128 0
	.uleb128 0x9
	.long	.LASF23
	.byte	0x8
	.value	0x115
	.long	0x120
	.sleb128 8
	.uleb128 0x9
	.long	.LASF24
	.byte	0x8
	.value	0x116
	.long	0x120
	.sleb128 16
	.uleb128 0x9
	.long	.LASF25
	.byte	0x8
	.value	0x117
	.long	0x120
	.sleb128 24
	.uleb128 0x9
	.long	.LASF26
	.byte	0x8
	.value	0x118
	.long	0x120
	.sleb128 32
	.uleb128 0x9
	.long	.LASF27
	.byte	0x8
	.value	0x119
	.long	0x120
	.sleb128 40
	.uleb128 0x9
	.long	.LASF28
	.byte	0x8
	.value	0x11a
	.long	0x120
	.sleb128 48
	.uleb128 0x9
	.long	.LASF29
	.byte	0x8
	.value	0x11b
	.long	0x120
	.sleb128 56
	.uleb128 0x9
	.long	.LASF30
	.byte	0x8
	.value	0x11c
	.long	0x120
	.sleb128 64
	.uleb128 0x9
	.long	.LASF31
	.byte	0x8
	.value	0x11e
	.long	0x120
	.sleb128 72
	.uleb128 0x9
	.long	.LASF32
	.byte	0x8
	.value	0x11f
	.long	0x120
	.sleb128 80
	.uleb128 0x9
	.long	.LASF33
	.byte	0x8
	.value	0x120
	.long	0x120
	.sleb128 88
	.uleb128 0x9
	.long	.LASF34
	.byte	0x8
	.value	0x122
	.long	0x301
	.sleb128 96
	.uleb128 0x9
	.long	.LASF35
	.byte	0x8
	.value	0x124
	.long	0x307
	.sleb128 104
	.uleb128 0x9
	.long	.LASF36
	.byte	0x8
	.value	0x126
	.long	0x2d
	.sleb128 112
	.uleb128 0x9
	.long	.LASF37
	.byte	0x8
	.value	0x12a
	.long	0x2d
	.sleb128 116
	.uleb128 0x9
	.long	.LASF38
	.byte	0x8
	.value	0x12c
	.long	0x10a
	.sleb128 120
	.uleb128 0x9
	.long	.LASF39
	.byte	0x8
	.value	0x130
	.long	0xee
	.sleb128 128
	.uleb128 0x9
	.long	.LASF40
	.byte	0x8
	.value	0x131
	.long	0xfc
	.sleb128 130
	.uleb128 0x9
	.long	.LASF41
	.byte	0x8
	.value	0x132
	.long	0x30d
	.sleb128 131
	.uleb128 0x9
	.long	.LASF42
	.byte	0x8
	.value	0x136
	.long	0x31d
	.sleb128 136
	.uleb128 0x9
	.long	.LASF43
	.byte	0x8
	.value	0x13f
	.long	0x115
	.sleb128 144
	.uleb128 0x9
	.long	.LASF44
	.byte	0x8
	.value	0x148
	.long	0x4b
	.sleb128 152
	.uleb128 0x9
	.long	.LASF45
	.byte	0x8
	.value	0x149
	.long	0x4b
	.sleb128 160
	.uleb128 0x9
	.long	.LASF46
	.byte	0x8
	.value	0x14a
	.long	0x4b
	.sleb128 168
	.uleb128 0x9
	.long	.LASF47
	.byte	0x8
	.value	0x14b
	.long	0x4b
	.sleb128 176
	.uleb128 0x9
	.long	.LASF48
	.byte	0x8
	.value	0x14c
	.long	0xd5
	.sleb128 184
	.uleb128 0x9
	.long	.LASF49
	.byte	0x8
	.value	0x14e
	.long	0x2d
	.sleb128 192
	.uleb128 0x9
	.long	.LASF50
	.byte	0x8
	.value	0x150
	.long	0x323
	.sleb128 196
	.byte	0x0
	.uleb128 0xa
	.long	.LASF751
	.byte	0x8
	.byte	0xb4
	.uleb128 0xb
	.long	.LASF51
	.byte	0x18
	.byte	0x8
	.byte	0xba
	.long	0x301
	.uleb128 0xc
	.long	.LASF52
	.byte	0x8
	.byte	0xbb
	.long	0x301
	.sleb128 0
	.uleb128 0xc
	.long	.LASF53
	.byte	0x8
	.byte	0xbc
	.long	0x307
	.sleb128 8
	.uleb128 0xc
	.long	.LASF54
	.byte	0x8
	.byte	0xc0
	.long	0x2d
	.sleb128 16
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2d0
	.uleb128 0x4
	.byte	0x8
	.long	0x12d
	.uleb128 0xd
	.long	0x126
	.long	0x31d
	.uleb128 0xe
	.long	0xe0
	.byte	0x0
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2c9
	.uleb128 0xd
	.long	0x126
	.long	0x333
	.uleb128 0xe
	.long	0xe0
	.byte	0x13
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x339
	.uleb128 0xf
	.long	0x126
	.uleb128 0x3
	.long	.LASF55
	.byte	0x2
	.byte	0x80
	.long	0x2d
	.uleb128 0x3
	.long	.LASF56
	.byte	0x2
	.byte	0x8d
	.long	0x2d
	.uleb128 0x3
	.long	.LASF57
	.byte	0x2
	.byte	0xab
	.long	0x2d
	.uleb128 0x3
	.long	.LASF58
	.byte	0x2
	.byte	0xbe
	.long	0x2d
	.uleb128 0x3
	.long	.LASF59
	.byte	0x2
	.byte	0xdc
	.long	0x2d
	.uleb128 0x10
	.byte	0x4
	.byte	0x2
	.byte	0xe9
	.long	0x38a
	.uleb128 0x11
	.long	.LASF60
	.sleb128 4
	.uleb128 0x11
	.long	.LASF61
	.sleb128 8
	.byte	0x0
	.uleb128 0x3
	.long	.LASF62
	.byte	0x2
	.byte	0xe9
	.long	0x375
	.uleb128 0x7
	.byte	0x4
	.byte	0x4
	.long	.LASF63
	.uleb128 0x3
	.long	.LASF64
	.byte	0x1
	.byte	0x23
	.long	0xce
	.uleb128 0x3
	.long	.LASF65
	.byte	0x1
	.byte	0x7e
	.long	0xce
	.uleb128 0x6
	.long	.LASF66
	.byte	0x1
	.value	0x146
	.long	0xce
	.uleb128 0x6
	.long	.LASF67
	.byte	0x1
	.value	0x151
	.long	0x3a7
	.uleb128 0x12
	.byte	0x4
	.byte	0x2
	.value	0x197
	.long	0x3e0
	.uleb128 0x11
	.long	.LASF68
	.sleb128 0
	.uleb128 0x11
	.long	.LASF69
	.sleb128 1
	.byte	0x0
	.uleb128 0x6
	.long	.LASF70
	.byte	0x2
	.value	0x197
	.long	0x3ca
	.uleb128 0x12
	.byte	0x4
	.byte	0x2
	.value	0x1a6
	.long	0x408
	.uleb128 0x11
	.long	.LASF71
	.sleb128 0
	.uleb128 0x11
	.long	.LASF72
	.sleb128 1
	.uleb128 0x11
	.long	.LASF73
	.sleb128 2
	.byte	0x0
	.uleb128 0x6
	.long	.LASF74
	.byte	0x2
	.value	0x4d7
	.long	0x414
	.uleb128 0x4
	.byte	0x8
	.long	0x41a
	.uleb128 0x13
	.long	.LASF75
	.value	0x1c0
	.byte	0x9
	.byte	0x35
	.long	0x6b4
	.uleb128 0xc
	.long	.LASF76
	.byte	0x9
	.byte	0x36
	.long	0x349
	.sleb128 0
	.uleb128 0xc
	.long	.LASF77
	.byte	0x9
	.byte	0x37
	.long	0x98b
	.sleb128 8
	.uleb128 0xc
	.long	.LASF78
	.byte	0x9
	.byte	0x38
	.long	0x34
	.sleb128 16
	.uleb128 0xc
	.long	.LASF79
	.byte	0x9
	.byte	0x39
	.long	0x36a
	.sleb128 20
	.uleb128 0xc
	.long	.LASF80
	.byte	0x9
	.byte	0x3a
	.long	0x3b2
	.sleb128 24
	.uleb128 0xc
	.long	.LASF81
	.byte	0x9
	.byte	0x3a
	.long	0x3b2
	.sleb128 32
	.uleb128 0x14
	.string	"mem"
	.byte	0x9
	.byte	0x3a
	.long	0x3b2
	.sleb128 40
	.uleb128 0x14
	.string	"id"
	.byte	0x9
	.byte	0x3b
	.long	0x36a
	.sleb128 48
	.uleb128 0xc
	.long	.LASF82
	.byte	0x9
	.byte	0x3c
	.long	0x36a
	.sleb128 52
	.uleb128 0x14
	.string	"tag"
	.byte	0x9
	.byte	0x3d
	.long	0x35f
	.sleb128 56
	.uleb128 0xc
	.long	.LASF83
	.byte	0x9
	.byte	0x3e
	.long	0x6b4
	.sleb128 64
	.uleb128 0xc
	.long	.LASF84
	.byte	0x9
	.byte	0x3f
	.long	0x80c
	.sleb128 72
	.uleb128 0xc
	.long	.LASF85
	.byte	0x9
	.byte	0x40
	.long	0x120
	.sleb128 80
	.uleb128 0xc
	.long	.LASF86
	.byte	0x9
	.byte	0x41
	.long	0x120
	.sleb128 88
	.uleb128 0xc
	.long	.LASF87
	.byte	0x9
	.byte	0x42
	.long	0x120
	.sleb128 96
	.uleb128 0xc
	.long	.LASF88
	.byte	0x9
	.byte	0x43
	.long	0x120
	.sleb128 104
	.uleb128 0xc
	.long	.LASF89
	.byte	0x9
	.byte	0x44
	.long	0x408
	.sleb128 112
	.uleb128 0xc
	.long	.LASF90
	.byte	0x9
	.byte	0x45
	.long	0x36a
	.sleb128 120
	.uleb128 0xc
	.long	.LASF91
	.byte	0x9
	.byte	0x46
	.long	0x120
	.sleb128 128
	.uleb128 0xc
	.long	.LASF92
	.byte	0x9
	.byte	0x47
	.long	0x120
	.sleb128 136
	.uleb128 0xc
	.long	.LASF93
	.byte	0x9
	.byte	0x48
	.long	0x36a
	.sleb128 144
	.uleb128 0x14
	.string	"cpp"
	.byte	0x9
	.byte	0x49
	.long	0x4b
	.sleb128 152
	.uleb128 0xc
	.long	.LASF94
	.byte	0x9
	.byte	0x4a
	.long	0x36a
	.sleb128 160
	.uleb128 0xc
	.long	.LASF95
	.byte	0x9
	.byte	0x4b
	.long	0x36a
	.sleb128 164
	.uleb128 0xc
	.long	.LASF96
	.byte	0x9
	.byte	0x4c
	.long	0x36a
	.sleb128 168
	.uleb128 0xc
	.long	.LASF97
	.byte	0x9
	.byte	0x4c
	.long	0x36a
	.sleb128 172
	.uleb128 0xc
	.long	.LASF98
	.byte	0x9
	.byte	0x4d
	.long	0x991
	.sleb128 176
	.uleb128 0xc
	.long	.LASF99
	.byte	0x9
	.byte	0x4d
	.long	0x991
	.sleb128 184
	.uleb128 0xc
	.long	.LASF100
	.byte	0x9
	.byte	0x4e
	.long	0x991
	.sleb128 192
	.uleb128 0xc
	.long	.LASF101
	.byte	0x9
	.byte	0x4e
	.long	0x997
	.sleb128 200
	.uleb128 0xc
	.long	.LASF102
	.byte	0x9
	.byte	0x4f
	.long	0x36a
	.sleb128 208
	.uleb128 0xc
	.long	.LASF103
	.byte	0x9
	.byte	0x4f
	.long	0x36a
	.sleb128 212
	.uleb128 0xc
	.long	.LASF104
	.byte	0x9
	.byte	0x50
	.long	0x991
	.sleb128 216
	.uleb128 0xc
	.long	.LASF105
	.byte	0x9
	.byte	0x50
	.long	0x991
	.sleb128 224
	.uleb128 0xc
	.long	.LASF106
	.byte	0x9
	.byte	0x51
	.long	0x824
	.sleb128 232
	.uleb128 0xc
	.long	.LASF107
	.byte	0x9
	.byte	0x51
	.long	0x99d
	.sleb128 240
	.uleb128 0xc
	.long	.LASF108
	.byte	0x9
	.byte	0x52
	.long	0x36a
	.sleb128 248
	.uleb128 0xc
	.long	.LASF109
	.byte	0x9
	.byte	0x52
	.long	0x36a
	.sleb128 252
	.uleb128 0xc
	.long	.LASF110
	.byte	0x9
	.byte	0x53
	.long	0x991
	.sleb128 256
	.uleb128 0xc
	.long	.LASF111
	.byte	0x9
	.byte	0x53
	.long	0x991
	.sleb128 264
	.uleb128 0xc
	.long	.LASF112
	.byte	0x9
	.byte	0x54
	.long	0x9a3
	.sleb128 272
	.uleb128 0xc
	.long	.LASF113
	.byte	0x9
	.byte	0x54
	.long	0x9a9
	.sleb128 280
	.uleb128 0xc
	.long	.LASF114
	.byte	0x9
	.byte	0x55
	.long	0x75a
	.sleb128 288
	.uleb128 0xc
	.long	.LASF115
	.byte	0x9
	.byte	0x56
	.long	0x4b
	.sleb128 296
	.uleb128 0xc
	.long	.LASF116
	.byte	0x9
	.byte	0x57
	.long	0x9bf
	.sleb128 304
	.uleb128 0xc
	.long	.LASF117
	.byte	0x9
	.byte	0x59
	.long	0x36a
	.sleb128 312
	.uleb128 0xc
	.long	.LASF118
	.byte	0x9
	.byte	0x5a
	.long	0x9c5
	.sleb128 320
	.uleb128 0xc
	.long	.LASF119
	.byte	0x9
	.byte	0x5b
	.long	0x9c5
	.sleb128 360
	.uleb128 0xc
	.long	.LASF120
	.byte	0x9
	.byte	0x5c
	.long	0x9f0
	.sleb128 400
	.uleb128 0xc
	.long	.LASF121
	.byte	0x9
	.byte	0x5d
	.long	0x38a
	.sleb128 440
	.uleb128 0xc
	.long	.LASF122
	.byte	0x9
	.byte	0x5e
	.long	0x3e0
	.sleb128 444
	.byte	0x0
	.uleb128 0x6
	.long	.LASF123
	.byte	0x2
	.value	0x4e1
	.long	0x6c0
	.uleb128 0x4
	.byte	0x8
	.long	0x6c6
	.uleb128 0x15
	.long	.LASF125
	.byte	0x1
	.uleb128 0x3
	.long	.LASF124
	.byte	0xa
	.byte	0x16
	.long	0x6d7
	.uleb128 0x4
	.byte	0x8
	.long	0x6dd
	.uleb128 0x15
	.long	.LASF126
	.byte	0x1
	.uleb128 0x3
	.long	.LASF127
	.byte	0xb
	.byte	0x99
	.long	0xf5
	.uleb128 0x7
	.byte	0x8
	.byte	0x7
	.long	.LASF128
	.uleb128 0xd
	.long	0x126
	.long	0x705
	.uleb128 0xe
	.long	0xe0
	.byte	0x1f
	.byte	0x0
	.uleb128 0x3
	.long	.LASF129
	.byte	0xc
	.byte	0x7b
	.long	0xe0
	.uleb128 0x4
	.byte	0x8
	.long	0x716
	.uleb128 0x16
	.uleb128 0x17
	.long	.LASF752
	.byte	0x4
	.byte	0x17
	.byte	0x75
	.long	0x75a
	.uleb128 0x11
	.long	.LASF130
	.sleb128 -1
	.uleb128 0x11
	.long	.LASF131
	.sleb128 0
	.uleb128 0x11
	.long	.LASF132
	.sleb128 1
	.uleb128 0x11
	.long	.LASF133
	.sleb128 2
	.uleb128 0x11
	.long	.LASF134
	.sleb128 3
	.uleb128 0x11
	.long	.LASF135
	.sleb128 4
	.uleb128 0x11
	.long	.LASF136
	.sleb128 5
	.uleb128 0x11
	.long	.LASF137
	.sleb128 6
	.uleb128 0x11
	.long	.LASF138
	.sleb128 7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x760
	.uleb128 0x4
	.byte	0x8
	.long	0x766
	.uleb128 0x18
	.byte	0x1
	.uleb128 0x4
	.byte	0x8
	.long	0x76e
	.uleb128 0x19
	.byte	0x1
	.long	0x33e
	.uleb128 0x12
	.byte	0x4
	.byte	0xd
	.value	0x164
	.long	0x790
	.uleb128 0x11
	.long	.LASF139
	.sleb128 0
	.uleb128 0x11
	.long	.LASF140
	.sleb128 1
	.uleb128 0x11
	.long	.LASF141
	.sleb128 2
	.byte	0x0
	.uleb128 0x1a
	.value	0x708
	.byte	0xd
	.value	0x185
	.long	0x7e0
	.uleb128 0x9
	.long	.LASF142
	.byte	0xd
	.value	0x186
	.long	0x7e0
	.sleb128 0
	.uleb128 0x9
	.long	.LASF143
	.byte	0xd
	.value	0x187
	.long	0x7e0
	.sleb128 512
	.uleb128 0x9
	.long	.LASF144
	.byte	0xd
	.value	0x188
	.long	0x7e0
	.sleb128 1024
	.uleb128 0x9
	.long	.LASF145
	.byte	0xd
	.value	0x189
	.long	0x7f0
	.sleb128 1536
	.uleb128 0x9
	.long	.LASF146
	.byte	0xd
	.value	0x18a
	.long	0x2d
	.sleb128 1792
	.byte	0x0
	.uleb128 0xd
	.long	0x333
	.long	0x7f0
	.uleb128 0xe
	.long	0xe0
	.byte	0x3f
	.byte	0x0
	.uleb128 0xd
	.long	0x2d
	.long	0x800
	.uleb128 0xe
	.long	0xe0
	.byte	0x3f
	.byte	0x0
	.uleb128 0x6
	.long	.LASF147
	.byte	0xd
	.value	0x18b
	.long	0x790
	.uleb128 0x6
	.long	.LASF148
	.byte	0x2
	.value	0x5b8
	.long	0x818
	.uleb128 0x4
	.byte	0x8
	.long	0x81e
	.uleb128 0x15
	.long	.LASF149
	.byte	0x1
	.uleb128 0x4
	.byte	0x8
	.long	0x39c
	.uleb128 0x1b
	.byte	0x40
	.byte	0x9
	.byte	0x23
	.long	0x893
	.uleb128 0xc
	.long	.LASF150
	.byte	0x9
	.byte	0x24
	.long	0x8a8
	.sleb128 0
	.uleb128 0xc
	.long	.LASF151
	.byte	0x9
	.byte	0x25
	.long	0x8c3
	.sleb128 8
	.uleb128 0xc
	.long	.LASF152
	.byte	0x9
	.byte	0x26
	.long	0x8df
	.sleb128 16
	.uleb128 0xc
	.long	.LASF153
	.byte	0x9
	.byte	0x27
	.long	0x8ff
	.sleb128 24
	.uleb128 0xc
	.long	.LASF154
	.byte	0x9
	.byte	0x28
	.long	0x91f
	.sleb128 32
	.uleb128 0xc
	.long	.LASF155
	.byte	0x9
	.byte	0x29
	.long	0x944
	.sleb128 40
	.uleb128 0xc
	.long	.LASF156
	.byte	0x9
	.byte	0x2a
	.long	0x964
	.sleb128 48
	.uleb128 0xc
	.long	.LASF157
	.byte	0x9
	.byte	0x2b
	.long	0x97a
	.sleb128 56
	.byte	0x0
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x8a8
	.uleb128 0x1d
	.long	0x408
	.uleb128 0x1d
	.long	0x3f
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x893
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x8c3
	.uleb128 0x1d
	.long	0x408
	.uleb128 0x1d
	.long	0x6cc
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x8ae
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x8d9
	.uleb128 0x1d
	.long	0x8d9
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x408
	.uleb128 0x4
	.byte	0x8
	.long	0x8c9
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x8ff
	.uleb128 0x1d
	.long	0x408
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x408
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x8e5
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x91f
	.uleb128 0x1d
	.long	0x408
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x8d9
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x905
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x944
	.uleb128 0x1d
	.long	0x408
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x760
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x925
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x964
	.uleb128 0x1d
	.long	0x408
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x75a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x94a
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x97a
	.uleb128 0x1d
	.long	0x408
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x96a
	.uleb128 0x3
	.long	.LASF158
	.byte	0x9
	.byte	0x2c
	.long	0x82a
	.uleb128 0x4
	.byte	0x8
	.long	0x980
	.uleb128 0x4
	.byte	0x8
	.long	0x36a
	.uleb128 0x4
	.byte	0x8
	.long	0x991
	.uleb128 0x4
	.byte	0x8
	.long	0x824
	.uleb128 0x4
	.byte	0x8
	.long	0x3a7
	.uleb128 0x4
	.byte	0x8
	.long	0x9a3
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x9bf
	.uleb128 0x1d
	.long	0x4b
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x9af
	.uleb128 0xd
	.long	0x9ea
	.long	0x9d5
	.uleb128 0xe
	.long	0xe0
	.byte	0x4
	.byte	0x0
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x9ea
	.uleb128 0x1d
	.long	0x408
	.uleb128 0x1d
	.long	0x4b
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x9d5
	.uleb128 0xd
	.long	0x4b
	.long	0xa00
	.uleb128 0xe
	.long	0xe0
	.byte	0x4
	.byte	0x0
	.uleb128 0x3
	.long	.LASF75
	.byte	0x9
	.byte	0x5f
	.long	0x41a
	.uleb128 0x3
	.long	.LASF159
	.byte	0xe
	.byte	0x43
	.long	0xa16
	.uleb128 0x4
	.byte	0x8
	.long	0xa1c
	.uleb128 0x15
	.long	.LASF160
	.byte	0x1
	.uleb128 0x1b
	.byte	0x10
	.byte	0xe
	.byte	0x4c
	.long	0xa43
	.uleb128 0xc
	.long	.LASF91
	.byte	0xe
	.byte	0x4d
	.long	0x120
	.sleb128 0
	.uleb128 0xc
	.long	.LASF76
	.byte	0xe
	.byte	0x4e
	.long	0x349
	.sleb128 8
	.byte	0x0
	.uleb128 0x3
	.long	.LASF161
	.byte	0xe
	.byte	0x4f
	.long	0xa22
	.uleb128 0x1b
	.byte	0x20
	.byte	0xe
	.byte	0x51
	.long	0xa92
	.uleb128 0x14
	.string	"id"
	.byte	0xe
	.byte	0x52
	.long	0x349
	.sleb128 0
	.uleb128 0xc
	.long	.LASF162
	.byte	0xe
	.byte	0x53
	.long	0x2d
	.sleb128 4
	.uleb128 0xc
	.long	.LASF163
	.byte	0xe
	.byte	0x54
	.long	0x2d
	.sleb128 8
	.uleb128 0x14
	.string	"mem"
	.byte	0xe
	.byte	0x55
	.long	0x3b2
	.sleb128 16
	.uleb128 0xc
	.long	.LASF164
	.byte	0xe
	.byte	0x56
	.long	0x3b2
	.sleb128 24
	.byte	0x0
	.uleb128 0x3
	.long	.LASF165
	.byte	0xe
	.byte	0x57
	.long	0xa4e
	.uleb128 0x3
	.long	.LASF166
	.byte	0xe
	.byte	0x59
	.long	0xaa8
	.uleb128 0x4
	.byte	0x8
	.long	0xaae
	.uleb128 0xb
	.long	.LASF167
	.byte	0x10
	.byte	0xe
	.byte	0x5a
	.long	0xadf
	.uleb128 0xc
	.long	.LASF168
	.byte	0xe
	.byte	0x5b
	.long	0x2d
	.sleb128 0
	.uleb128 0xc
	.long	.LASF169
	.byte	0xe
	.byte	0x5c
	.long	0x2d
	.sleb128 4
	.uleb128 0xc
	.long	.LASF170
	.byte	0xe
	.byte	0x5d
	.long	0xadf
	.sleb128 8
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xa43
	.uleb128 0x3
	.long	.LASF171
	.byte	0xe
	.byte	0x60
	.long	0xaf0
	.uleb128 0x4
	.byte	0x8
	.long	0xaf6
	.uleb128 0xb
	.long	.LASF172
	.byte	0x10
	.byte	0xe
	.byte	0x61
	.long	0xb27
	.uleb128 0xc
	.long	.LASF168
	.byte	0xe
	.byte	0x62
	.long	0x2d
	.sleb128 0
	.uleb128 0xc
	.long	.LASF169
	.byte	0xe
	.byte	0x63
	.long	0x2d
	.sleb128 4
	.uleb128 0xc
	.long	.LASF170
	.byte	0xe
	.byte	0x64
	.long	0xb27
	.sleb128 8
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xa92
	.uleb128 0x1b
	.byte	0x10
	.byte	0xe
	.byte	0x70
	.long	0xb4e
	.uleb128 0xc
	.long	.LASF91
	.byte	0xe
	.byte	0x71
	.long	0x120
	.sleb128 0
	.uleb128 0xc
	.long	.LASF76
	.byte	0xe
	.byte	0x72
	.long	0x349
	.sleb128 8
	.byte	0x0
	.uleb128 0x3
	.long	.LASF173
	.byte	0xe
	.byte	0x77
	.long	0xb2d
	.uleb128 0x1b
	.byte	0x40
	.byte	0xe
	.byte	0x79
	.long	0xbd9
	.uleb128 0x14
	.string	"id"
	.byte	0xe
	.byte	0x7a
	.long	0x2d
	.sleb128 0
	.uleb128 0xc
	.long	.LASF174
	.byte	0xe
	.byte	0x7b
	.long	0x3e0
	.sleb128 4
	.uleb128 0xc
	.long	.LASF175
	.byte	0xe
	.byte	0x7c
	.long	0x3e0
	.sleb128 8
	.uleb128 0xc
	.long	.LASF176
	.byte	0xe
	.byte	0x7d
	.long	0x2d
	.sleb128 12
	.uleb128 0xc
	.long	.LASF4
	.byte	0xe
	.byte	0x7e
	.long	0x2d
	.sleb128 16
	.uleb128 0xc
	.long	.LASF80
	.byte	0xe
	.byte	0x7f
	.long	0x3b2
	.sleb128 24
	.uleb128 0xc
	.long	.LASF81
	.byte	0xe
	.byte	0x80
	.long	0x3b2
	.sleb128 32
	.uleb128 0xc
	.long	.LASF177
	.byte	0xe
	.byte	0x81
	.long	0x3b2
	.sleb128 40
	.uleb128 0xc
	.long	.LASF178
	.byte	0xe
	.byte	0x82
	.long	0x3b2
	.sleb128 48
	.uleb128 0xc
	.long	.LASF179
	.byte	0xe
	.byte	0x83
	.long	0x3b2
	.sleb128 56
	.byte	0x0
	.uleb128 0x3
	.long	.LASF180
	.byte	0xe
	.byte	0x84
	.long	0xb59
	.uleb128 0x3
	.long	.LASF181
	.byte	0xe
	.byte	0x86
	.long	0xbef
	.uleb128 0x4
	.byte	0x8
	.long	0xbf5
	.uleb128 0xb
	.long	.LASF182
	.byte	0x10
	.byte	0xe
	.byte	0x87
	.long	0xc26
	.uleb128 0xc
	.long	.LASF183
	.byte	0xe
	.byte	0x88
	.long	0x2d
	.sleb128 0
	.uleb128 0xc
	.long	.LASF184
	.byte	0xe
	.byte	0x89
	.long	0x2d
	.sleb128 4
	.uleb128 0xc
	.long	.LASF185
	.byte	0xe
	.byte	0x8a
	.long	0xc26
	.sleb128 8
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xb4e
	.uleb128 0x3
	.long	.LASF186
	.byte	0xe
	.byte	0x8d
	.long	0xc37
	.uleb128 0x4
	.byte	0x8
	.long	0xc3d
	.uleb128 0xb
	.long	.LASF187
	.byte	0x10
	.byte	0xe
	.byte	0x8e
	.long	0xc6e
	.uleb128 0xc
	.long	.LASF183
	.byte	0xe
	.byte	0x8f
	.long	0x2d
	.sleb128 0
	.uleb128 0xc
	.long	.LASF184
	.byte	0xe
	.byte	0x90
	.long	0x2d
	.sleb128 4
	.uleb128 0xc
	.long	.LASF185
	.byte	0xe
	.byte	0x91
	.long	0xc6e
	.sleb128 8
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xbd9
	.uleb128 0xb
	.long	.LASF188
	.byte	0x60
	.byte	0xe
	.byte	0x99
	.long	0xcbf
	.uleb128 0xc
	.long	.LASF91
	.byte	0xe
	.byte	0x9a
	.long	0x120
	.sleb128 0
	.uleb128 0xc
	.long	.LASF189
	.byte	0xe
	.byte	0x9b
	.long	0x3e0
	.sleb128 8
	.uleb128 0xc
	.long	.LASF190
	.byte	0xe
	.byte	0x9c
	.long	0xbd9
	.sleb128 16
	.uleb128 0xc
	.long	.LASF191
	.byte	0xe
	.byte	0x9d
	.long	0xc2c
	.sleb128 80
	.uleb128 0xc
	.long	.LASF192
	.byte	0xe
	.byte	0x9e
	.long	0xae5
	.sleb128 88
	.byte	0x0
	.uleb128 0x3
	.long	.LASF193
	.byte	0xe
	.byte	0x9f
	.long	0xc74
	.uleb128 0x3
	.long	.LASF194
	.byte	0xe
	.byte	0xa1
	.long	0xcd5
	.uleb128 0x4
	.byte	0x8
	.long	0xcdb
	.uleb128 0xb
	.long	.LASF195
	.byte	0x30
	.byte	0xe
	.byte	0xa3
	.long	0xd3c
	.uleb128 0xc
	.long	.LASF196
	.byte	0xe
	.byte	0xa4
	.long	0x2d
	.sleb128 0
	.uleb128 0xc
	.long	.LASF197
	.byte	0xe
	.byte	0xa5
	.long	0x2d
	.sleb128 4
	.uleb128 0xc
	.long	.LASF198
	.byte	0xe
	.byte	0xa6
	.long	0xa0b
	.sleb128 8
	.uleb128 0xc
	.long	.LASF199
	.byte	0xe
	.byte	0xa7
	.long	0x2d
	.sleb128 16
	.uleb128 0xc
	.long	.LASF200
	.byte	0xe
	.byte	0xa8
	.long	0xd3c
	.sleb128 24
	.uleb128 0xc
	.long	.LASF191
	.byte	0xe
	.byte	0xa9
	.long	0xbe4
	.sleb128 32
	.uleb128 0xc
	.long	.LASF192
	.byte	0xe
	.byte	0xaa
	.long	0xa9d
	.sleb128 40
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xcbf
	.uleb128 0x6
	.long	.LASF201
	.byte	0x2
	.value	0x8d2
	.long	0xd4e
	.uleb128 0x4
	.byte	0x8
	.long	0xd54
	.uleb128 0x15
	.long	.LASF202
	.byte	0x1
	.uleb128 0x12
	.byte	0x4
	.byte	0x2
	.value	0x963
	.long	0xd88
	.uleb128 0x11
	.long	.LASF203
	.sleb128 0
	.uleb128 0x11
	.long	.LASF204
	.sleb128 1
	.uleb128 0x11
	.long	.LASF205
	.sleb128 2
	.uleb128 0x11
	.long	.LASF206
	.sleb128 3
	.uleb128 0x11
	.long	.LASF207
	.sleb128 4
	.uleb128 0x11
	.long	.LASF208
	.sleb128 5
	.byte	0x0
	.uleb128 0x6
	.long	.LASF209
	.byte	0x2
	.value	0x963
	.long	0xd5a
	.uleb128 0x1e
	.string	"IS"
	.byte	0xf
	.byte	0x18
	.long	0xd9e
	.uleb128 0x4
	.byte	0x8
	.long	0xda4
	.uleb128 0x15
	.long	.LASF210
	.byte	0x1
	.uleb128 0x13
	.long	.LASF211
	.value	0x1e8
	.byte	0xf
	.byte	0xae
	.long	0xe10
	.uleb128 0x14
	.string	"hdr"
	.byte	0xf
	.byte	0xaf
	.long	0xa00
	.sleb128 0
	.uleb128 0x14
	.string	"ops"
	.byte	0xf
	.byte	0xaf
	.long	0x45
	.sleb128 448
	.uleb128 0x14
	.string	"n"
	.byte	0xf
	.byte	0xb0
	.long	0x36a
	.sleb128 456
	.uleb128 0xc
	.long	.LASF212
	.byte	0xf
	.byte	0xb1
	.long	0x991
	.sleb128 464
	.uleb128 0xc
	.long	.LASF213
	.byte	0xf
	.byte	0xb2
	.long	0x36a
	.sleb128 472
	.uleb128 0xc
	.long	.LASF214
	.byte	0xf
	.byte	0xb3
	.long	0x36a
	.sleb128 476
	.uleb128 0xc
	.long	.LASF215
	.byte	0xf
	.byte	0xb4
	.long	0x991
	.sleb128 480
	.byte	0x0
	.uleb128 0x3
	.long	.LASF216
	.byte	0xf
	.byte	0xb6
	.long	0xe1b
	.uleb128 0x4
	.byte	0x8
	.long	0xdaa
	.uleb128 0x10
	.byte	0x4
	.byte	0xf
	.byte	0xf4
	.long	0xe36
	.uleb128 0x11
	.long	.LASF217
	.sleb128 0
	.uleb128 0x11
	.long	.LASF218
	.sleb128 1
	.byte	0x0
	.uleb128 0x3
	.long	.LASF219
	.byte	0xf
	.byte	0xf4
	.long	0xe21
	.uleb128 0x3
	.long	.LASF220
	.byte	0xf
	.byte	0xf6
	.long	0xee
	.uleb128 0x8
	.long	.LASF221
	.byte	0x28
	.byte	0xf
	.value	0x106
	.long	0xeb0
	.uleb128 0x9
	.long	.LASF82
	.byte	0xf
	.value	0x107
	.long	0x36a
	.sleb128 0
	.uleb128 0x1f
	.string	"n"
	.byte	0xf
	.value	0x108
	.long	0x36a
	.sleb128 4
	.uleb128 0x1f
	.string	"is"
	.byte	0xf
	.value	0x109
	.long	0xeb0
	.sleb128 8
	.uleb128 0x9
	.long	.LASF78
	.byte	0xf
	.value	0x10a
	.long	0x34
	.sleb128 16
	.uleb128 0x9
	.long	.LASF222
	.byte	0xf
	.value	0x10b
	.long	0xeb6
	.sleb128 24
	.uleb128 0x1f
	.string	"N"
	.byte	0xf
	.value	0x10c
	.long	0x36a
	.sleb128 32
	.uleb128 0x9
	.long	.LASF223
	.byte	0xf
	.value	0x10d
	.long	0xe36
	.sleb128 36
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xd94
	.uleb128 0x4
	.byte	0x8
	.long	0xe41
	.uleb128 0x6
	.long	.LASF224
	.byte	0xf
	.value	0x10f
	.long	0xec8
	.uleb128 0x4
	.byte	0x8
	.long	0xe4c
	.uleb128 0x1e
	.string	"Vec"
	.byte	0x10
	.byte	0x16
	.long	0xed9
	.uleb128 0x4
	.byte	0x8
	.long	0xedf
	.uleb128 0x13
	.long	.LASF225
	.value	0x340
	.byte	0x3
	.byte	0xfa
	.long	0xf71
	.uleb128 0x14
	.string	"hdr"
	.byte	0x3
	.byte	0xfb
	.long	0xa00
	.sleb128 0
	.uleb128 0x14
	.string	"ops"
	.byte	0x3
	.byte	0xfb
	.long	0x1056
	.sleb128 448
	.uleb128 0x14
	.string	"map"
	.byte	0x3
	.byte	0xfc
	.long	0xfc5
	.sleb128 456
	.uleb128 0xc
	.long	.LASF226
	.byte	0x3
	.byte	0xfd
	.long	0x4b
	.sleb128 464
	.uleb128 0xc
	.long	.LASF227
	.byte	0x3
	.byte	0xfe
	.long	0x3e0
	.sleb128 472
	.uleb128 0xc
	.long	.LASF228
	.byte	0x3
	.byte	0xff
	.long	0x18e7
	.sleb128 480
	.uleb128 0xc
	.long	.LASF229
	.byte	0x3
	.byte	0xff
	.long	0x18e7
	.sleb128 648
	.uleb128 0x9
	.long	.LASF230
	.byte	0x3
	.value	0x100
	.long	0x3e0
	.sleb128 816
	.uleb128 0x9
	.long	.LASF231
	.byte	0x3
	.value	0x102
	.long	0x1913
	.sleb128 820
	.uleb128 0x9
	.long	.LASF232
	.byte	0x3
	.value	0x103
	.long	0x4b
	.sleb128 824
	.byte	0x0
	.uleb128 0x10
	.byte	0x4
	.byte	0x10
	.byte	0x9f
	.long	0xf98
	.uleb128 0x11
	.long	.LASF233
	.sleb128 0
	.uleb128 0x11
	.long	.LASF234
	.sleb128 1
	.uleb128 0x11
	.long	.LASF235
	.sleb128 2
	.uleb128 0x11
	.long	.LASF236
	.sleb128 3
	.uleb128 0x11
	.long	.LASF237
	.sleb128 4
	.byte	0x0
	.uleb128 0x3
	.long	.LASF238
	.byte	0x10
	.byte	0x9f
	.long	0xf71
	.uleb128 0x12
	.byte	0x4
	.byte	0x10
	.value	0x1e7
	.long	0xfb9
	.uleb128 0x11
	.long	.LASF239
	.sleb128 0
	.uleb128 0x11
	.long	.LASF240
	.sleb128 1
	.byte	0x0
	.uleb128 0x6
	.long	.LASF241
	.byte	0x10
	.value	0x1e7
	.long	0xfa3
	.uleb128 0x3
	.long	.LASF242
	.byte	0x3
	.byte	0x16
	.long	0xfd0
	.uleb128 0x4
	.byte	0x8
	.long	0xfd6
	.uleb128 0xb
	.long	.LASF243
	.byte	0x38
	.byte	0x3
	.byte	0x17
	.long	0x1056
	.uleb128 0xc
	.long	.LASF78
	.byte	0x3
	.byte	0x18
	.long	0x34
	.sleb128 0
	.uleb128 0x14
	.string	"n"
	.byte	0x3
	.byte	0x19
	.long	0x36a
	.sleb128 4
	.uleb128 0x14
	.string	"N"
	.byte	0x3
	.byte	0x19
	.long	0x36a
	.sleb128 8
	.uleb128 0xc
	.long	.LASF244
	.byte	0x3
	.byte	0x1a
	.long	0x36a
	.sleb128 12
	.uleb128 0xc
	.long	.LASF245
	.byte	0x3
	.byte	0x1a
	.long	0x36a
	.sleb128 16
	.uleb128 0xc
	.long	.LASF246
	.byte	0x3
	.byte	0x1b
	.long	0x991
	.sleb128 24
	.uleb128 0x14
	.string	"bs"
	.byte	0x3
	.byte	0x1c
	.long	0x36a
	.sleb128 32
	.uleb128 0xc
	.long	.LASF247
	.byte	0x3
	.byte	0x1d
	.long	0x36a
	.sleb128 36
	.uleb128 0xc
	.long	.LASF248
	.byte	0x3
	.byte	0x1e
	.long	0xe10
	.sleb128 40
	.uleb128 0xc
	.long	.LASF249
	.byte	0x3
	.byte	0x1f
	.long	0xe10
	.sleb128 48
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x105c
	.uleb128 0x13
	.long	.LASF250
	.value	0x200
	.byte	0x3
	.byte	0x8f
	.long	0x13a2
	.uleb128 0xc
	.long	.LASF251
	.byte	0x3
	.byte	0x90
	.long	0x13bd
	.sleb128 0
	.uleb128 0xc
	.long	.LASF252
	.byte	0x3
	.byte	0x91
	.long	0x13e3
	.sleb128 8
	.uleb128 0xc
	.long	.LASF253
	.byte	0x3
	.byte	0x92
	.long	0x13fe
	.sleb128 16
	.uleb128 0x14
	.string	"dot"
	.byte	0x3
	.byte	0x93
	.long	0x141e
	.sleb128 24
	.uleb128 0xc
	.long	.LASF254
	.byte	0x3
	.byte	0x94
	.long	0x144e
	.sleb128 32
	.uleb128 0xc
	.long	.LASF255
	.byte	0x3
	.byte	0x95
	.long	0x146e
	.sleb128 40
	.uleb128 0xc
	.long	.LASF256
	.byte	0x3
	.byte	0x96
	.long	0x141e
	.sleb128 48
	.uleb128 0xc
	.long	.LASF257
	.byte	0x3
	.byte	0x97
	.long	0x144e
	.sleb128 56
	.uleb128 0xc
	.long	.LASF258
	.byte	0x3
	.byte	0x98
	.long	0x1489
	.sleb128 64
	.uleb128 0xc
	.long	.LASF259
	.byte	0x3
	.byte	0x99
	.long	0x14a4
	.sleb128 72
	.uleb128 0x14
	.string	"set"
	.byte	0x3
	.byte	0x9a
	.long	0x1489
	.sleb128 80
	.uleb128 0xc
	.long	.LASF260
	.byte	0x3
	.byte	0x9b
	.long	0x14a4
	.sleb128 88
	.uleb128 0xc
	.long	.LASF261
	.byte	0x3
	.byte	0x9c
	.long	0x14c4
	.sleb128 96
	.uleb128 0xc
	.long	.LASF262
	.byte	0x3
	.byte	0x9d
	.long	0x14e9
	.sleb128 104
	.uleb128 0xc
	.long	.LASF263
	.byte	0x3
	.byte	0x9e
	.long	0x1519
	.sleb128 112
	.uleb128 0xc
	.long	.LASF264
	.byte	0x3
	.byte	0x9f
	.long	0x14c4
	.sleb128 120
	.uleb128 0xc
	.long	.LASF265
	.byte	0x3
	.byte	0xa0
	.long	0x153e
	.sleb128 128
	.uleb128 0xc
	.long	.LASF266
	.byte	0x3
	.byte	0xa1
	.long	0x156d
	.sleb128 136
	.uleb128 0xc
	.long	.LASF267
	.byte	0x3
	.byte	0xa2
	.long	0x158d
	.sleb128 144
	.uleb128 0xc
	.long	.LASF268
	.byte	0x3
	.byte	0xa3
	.long	0x158d
	.sleb128 152
	.uleb128 0xc
	.long	.LASF269
	.byte	0x3
	.byte	0xa4
	.long	0x15c2
	.sleb128 160
	.uleb128 0xc
	.long	.LASF270
	.byte	0x3
	.byte	0xa5
	.long	0x15d8
	.sleb128 168
	.uleb128 0xc
	.long	.LASF271
	.byte	0x3
	.byte	0xa6
	.long	0x15d8
	.sleb128 176
	.uleb128 0xc
	.long	.LASF272
	.byte	0x3
	.byte	0xa7
	.long	0x15f3
	.sleb128 184
	.uleb128 0xc
	.long	.LASF273
	.byte	0x3
	.byte	0xa8
	.long	0x160e
	.sleb128 192
	.uleb128 0xc
	.long	.LASF274
	.byte	0x3
	.byte	0xa9
	.long	0x160e
	.sleb128 200
	.uleb128 0xc
	.long	.LASF275
	.byte	0x3
	.byte	0xaa
	.long	0x15f3
	.sleb128 208
	.uleb128 0x14
	.string	"max"
	.byte	0x3
	.byte	0xab
	.long	0x162e
	.sleb128 216
	.uleb128 0x14
	.string	"min"
	.byte	0x3
	.byte	0xac
	.long	0x162e
	.sleb128 224
	.uleb128 0xc
	.long	.LASF276
	.byte	0x3
	.byte	0xad
	.long	0x1649
	.sleb128 232
	.uleb128 0xc
	.long	.LASF277
	.byte	0x3
	.byte	0xae
	.long	0x1669
	.sleb128 240
	.uleb128 0xc
	.long	.LASF278
	.byte	0x3
	.byte	0xaf
	.long	0x15c2
	.sleb128 248
	.uleb128 0xc
	.long	.LASF152
	.byte	0x3
	.byte	0xb0
	.long	0x15d8
	.sleb128 256
	.uleb128 0xc
	.long	.LASF151
	.byte	0x3
	.byte	0xb1
	.long	0x1684
	.sleb128 264
	.uleb128 0xc
	.long	.LASF279
	.byte	0x3
	.byte	0xb2
	.long	0x169f
	.sleb128 272
	.uleb128 0xc
	.long	.LASF280
	.byte	0x3
	.byte	0xb3
	.long	0x169f
	.sleb128 280
	.uleb128 0xc
	.long	.LASF281
	.byte	0x3
	.byte	0xb4
	.long	0x141e
	.sleb128 288
	.uleb128 0xc
	.long	.LASF282
	.byte	0x3
	.byte	0xb5
	.long	0x141e
	.sleb128 296
	.uleb128 0xc
	.long	.LASF283
	.byte	0x3
	.byte	0xb6
	.long	0x146e
	.sleb128 304
	.uleb128 0xc
	.long	.LASF284
	.byte	0x3
	.byte	0xb7
	.long	0x144e
	.sleb128 312
	.uleb128 0xc
	.long	.LASF285
	.byte	0x3
	.byte	0xb8
	.long	0x144e
	.sleb128 320
	.uleb128 0xc
	.long	.LASF286
	.byte	0x3
	.byte	0xb9
	.long	0x1684
	.sleb128 328
	.uleb128 0xc
	.long	.LASF287
	.byte	0x3
	.byte	0xba
	.long	0x15d8
	.sleb128 336
	.uleb128 0xc
	.long	.LASF288
	.byte	0x3
	.byte	0xbb
	.long	0x15d8
	.sleb128 344
	.uleb128 0xc
	.long	.LASF289
	.byte	0x3
	.byte	0xbc
	.long	0x16ba
	.sleb128 352
	.uleb128 0xc
	.long	.LASF290
	.byte	0x3
	.byte	0xbd
	.long	0x15c2
	.sleb128 360
	.uleb128 0xc
	.long	.LASF291
	.byte	0x3
	.byte	0xbe
	.long	0x15d8
	.sleb128 368
	.uleb128 0xc
	.long	.LASF292
	.byte	0x3
	.byte	0xbf
	.long	0x15d8
	.sleb128 376
	.uleb128 0xc
	.long	.LASF293
	.byte	0x3
	.byte	0xc0
	.long	0x16da
	.sleb128 384
	.uleb128 0xc
	.long	.LASF294
	.byte	0x3
	.byte	0xc1
	.long	0x158d
	.sleb128 392
	.uleb128 0xc
	.long	.LASF295
	.byte	0x3
	.byte	0xc2
	.long	0x158d
	.sleb128 400
	.uleb128 0xc
	.long	.LASF296
	.byte	0x3
	.byte	0xc3
	.long	0x158d
	.sleb128 408
	.uleb128 0xc
	.long	.LASF297
	.byte	0x3
	.byte	0xc4
	.long	0x16ff
	.sleb128 416
	.uleb128 0xc
	.long	.LASF298
	.byte	0x3
	.byte	0xc5
	.long	0x15d8
	.sleb128 424
	.uleb128 0x14
	.string	"abs"
	.byte	0x3
	.byte	0xc6
	.long	0x15d8
	.sleb128 432
	.uleb128 0x14
	.string	"exp"
	.byte	0x3
	.byte	0xc7
	.long	0x15d8
	.sleb128 440
	.uleb128 0x14
	.string	"log"
	.byte	0x3
	.byte	0xc8
	.long	0x15d8
	.sleb128 448
	.uleb128 0xc
	.long	.LASF299
	.byte	0x3
	.byte	0xc9
	.long	0x15d8
	.sleb128 456
	.uleb128 0xc
	.long	.LASF300
	.byte	0x3
	.byte	0xca
	.long	0x15d8
	.sleb128 464
	.uleb128 0xc
	.long	.LASF301
	.byte	0x3
	.byte	0xcb
	.long	0x1724
	.sleb128 472
	.uleb128 0xc
	.long	.LASF302
	.byte	0x3
	.byte	0xcc
	.long	0x1724
	.sleb128 480
	.uleb128 0xc
	.long	.LASF303
	.byte	0x3
	.byte	0xcd
	.long	0x1749
	.sleb128 488
	.uleb128 0xc
	.long	.LASF304
	.byte	0x3
	.byte	0xce
	.long	0x1769
	.sleb128 496
	.uleb128 0xc
	.long	.LASF305
	.byte	0x3
	.byte	0xcf
	.long	0x1769
	.sleb128 504
	.byte	0x0
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x13b7
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x13b7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xece
	.uleb128 0x4
	.byte	0x8
	.long	0x13a2
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x13dd
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x13dd
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x13b7
	.uleb128 0x4
	.byte	0x8
	.long	0x13c3
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x13fe
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x13b7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x13e9
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x141e
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x9a3
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1404
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1443
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x1443
	.uleb128 0x1d
	.long	0x9a3
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1449
	.uleb128 0xf
	.long	0xece
	.uleb128 0x4
	.byte	0x8
	.long	0x1424
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x146e
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xf98
	.uleb128 0x1d
	.long	0x824
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1454
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1489
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x3a7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1474
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x14a4
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x148f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x14c4
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x14aa
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x14e9
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x14ca
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x150e
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x150e
	.uleb128 0x1d
	.long	0x13b7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1514
	.uleb128 0xf
	.long	0x3a7
	.uleb128 0x4
	.byte	0x8
	.long	0x14ef
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x153e
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x151f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x156d
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1544
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x158d
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1573
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x15b7
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x150e
	.uleb128 0x1d
	.long	0xd88
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x15bd
	.uleb128 0xf
	.long	0x36a
	.uleb128 0x4
	.byte	0x8
	.long	0x1593
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x15d8
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x15c8
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x15f3
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x9a9
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x15de
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x160e
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x991
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x15f9
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x162e
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x991
	.uleb128 0x1d
	.long	0x824
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1614
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1649
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xd42
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1634
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1669
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xfb9
	.uleb128 0x1d
	.long	0x3e0
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x164f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1684
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x6cc
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x166f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x169f
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x150e
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x168a
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x16ba
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xe10
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x16a5
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x16da
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x824
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x16c0
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x16ff
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x9a3
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x16e0
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1724
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xd88
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1705
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1749
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x9a3
	.uleb128 0x1d
	.long	0x9a3
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x172a
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x1769
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0x13b7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x174f
	.uleb128 0x1b
	.byte	0xa8
	.byte	0x3
	.byte	0xd8
	.long	0x18e1
	.uleb128 0xc
	.long	.LASF306
	.byte	0x3
	.byte	0xd9
	.long	0x36a
	.sleb128 0
	.uleb128 0xc
	.long	.LASF307
	.byte	0x3
	.byte	0xda
	.long	0x36a
	.sleb128 4
	.uleb128 0xc
	.long	.LASF308
	.byte	0x3
	.byte	0xdb
	.long	0x36a
	.sleb128 8
	.uleb128 0x14
	.string	"n"
	.byte	0x3
	.byte	0xdc
	.long	0x36a
	.sleb128 12
	.uleb128 0x14
	.string	"bs"
	.byte	0x3
	.byte	0xdd
	.long	0x36a
	.sleb128 16
	.uleb128 0xc
	.long	.LASF309
	.byte	0x3
	.byte	0xde
	.long	0x36a
	.sleb128 20
	.uleb128 0x14
	.string	"idx"
	.byte	0x3
	.byte	0xdf
	.long	0x991
	.sleb128 24
	.uleb128 0xc
	.long	.LASF310
	.byte	0x3
	.byte	0xe0
	.long	0x9a3
	.sleb128 32
	.uleb128 0xc
	.long	.LASF78
	.byte	0x3
	.byte	0xe2
	.long	0x34
	.sleb128 40
	.uleb128 0xc
	.long	.LASF311
	.byte	0x3
	.byte	0xe3
	.long	0x35f
	.sleb128 44
	.uleb128 0xc
	.long	.LASF312
	.byte	0x3
	.byte	0xe3
	.long	0x35f
	.sleb128 48
	.uleb128 0xc
	.long	.LASF313
	.byte	0x3
	.byte	0xe4
	.long	0x35f
	.sleb128 52
	.uleb128 0xc
	.long	.LASF314
	.byte	0x3
	.byte	0xe4
	.long	0x35f
	.sleb128 56
	.uleb128 0xc
	.long	.LASF315
	.byte	0x3
	.byte	0xe5
	.long	0x18e1
	.sleb128 64
	.uleb128 0xc
	.long	.LASF316
	.byte	0x3
	.byte	0xe6
	.long	0x18e1
	.sleb128 72
	.uleb128 0xc
	.long	.LASF317
	.byte	0x3
	.byte	0xe7
	.long	0xc2
	.sleb128 80
	.uleb128 0xc
	.long	.LASF318
	.byte	0x3
	.byte	0xe8
	.long	0x36a
	.sleb128 88
	.uleb128 0xc
	.long	.LASF319
	.byte	0x3
	.byte	0xe8
	.long	0x36a
	.sleb128 92
	.uleb128 0xc
	.long	.LASF320
	.byte	0x3
	.byte	0xe9
	.long	0x9a3
	.sleb128 96
	.uleb128 0xc
	.long	.LASF321
	.byte	0x3
	.byte	0xe9
	.long	0x9a3
	.sleb128 104
	.uleb128 0xc
	.long	.LASF322
	.byte	0x3
	.byte	0xea
	.long	0x991
	.sleb128 112
	.uleb128 0xc
	.long	.LASF323
	.byte	0x3
	.byte	0xea
	.long	0x991
	.sleb128 120
	.uleb128 0xc
	.long	.LASF324
	.byte	0x3
	.byte	0xeb
	.long	0x36a
	.sleb128 128
	.uleb128 0xc
	.long	.LASF325
	.byte	0x3
	.byte	0xec
	.long	0x991
	.sleb128 136
	.uleb128 0xc
	.long	.LASF326
	.byte	0x3
	.byte	0xed
	.long	0x36a
	.sleb128 144
	.uleb128 0xc
	.long	.LASF327
	.byte	0x3
	.byte	0xee
	.long	0x3e0
	.sleb128 148
	.uleb128 0xc
	.long	.LASF328
	.byte	0x3
	.byte	0xef
	.long	0x3e0
	.sleb128 152
	.uleb128 0xc
	.long	.LASF329
	.byte	0x3
	.byte	0xf0
	.long	0xd88
	.sleb128 156
	.uleb128 0xc
	.long	.LASF330
	.byte	0x3
	.byte	0xf1
	.long	0x991
	.sleb128 160
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x4d
	.uleb128 0x3
	.long	.LASF331
	.byte	0x3
	.byte	0xf2
	.long	0x176f
	.uleb128 0x10
	.byte	0x4
	.byte	0x3
	.byte	0xf5
	.long	0x1913
	.uleb128 0x11
	.long	.LASF332
	.sleb128 0
	.uleb128 0x11
	.long	.LASF333
	.sleb128 1
	.uleb128 0x11
	.long	.LASF334
	.sleb128 2
	.uleb128 0x11
	.long	.LASF335
	.sleb128 3
	.byte	0x0
	.uleb128 0x3
	.long	.LASF336
	.byte	0x3
	.byte	0xf5
	.long	0x18f2
	.uleb128 0x4
	.byte	0x8
	.long	0x35f
	.uleb128 0x8
	.long	.LASF337
	.byte	0x10
	.byte	0x10
	.value	0x21f
	.long	0x1948
	.uleb128 0x1f
	.string	"n"
	.byte	0x10
	.value	0x21f
	.long	0x36a
	.sleb128 0
	.uleb128 0x1f
	.string	"v"
	.byte	0x10
	.value	0x21f
	.long	0xece
	.sleb128 8
	.byte	0x0
	.uleb128 0x6
	.long	.LASF338
	.byte	0x10
	.value	0x220
	.long	0x1954
	.uleb128 0x4
	.byte	0x8
	.long	0x1924
	.uleb128 0x1e
	.string	"Mat"
	.byte	0x11
	.byte	0x12
	.long	0x1965
	.uleb128 0x4
	.byte	0x8
	.long	0x196b
	.uleb128 0x20
	.long	.LASF339
	.value	0x400
	.byte	0x12
	.value	0x11b
	.long	0x1b39
	.uleb128 0x1f
	.string	"hdr"
	.byte	0x12
	.value	0x11c
	.long	0xa00
	.sleb128 0
	.uleb128 0x1f
	.string	"ops"
	.byte	0x12
	.value	0x11c
	.long	0x2093
	.sleb128 448
	.uleb128 0x9
	.long	.LASF340
	.byte	0x12
	.value	0x11d
	.long	0xfc5
	.sleb128 456
	.uleb128 0x9
	.long	.LASF341
	.byte	0x12
	.value	0x11d
	.long	0xfc5
	.sleb128 464
	.uleb128 0x9
	.long	.LASF226
	.byte	0x12
	.value	0x11e
	.long	0x4b
	.sleb128 472
	.uleb128 0x9
	.long	.LASF342
	.byte	0x12
	.value	0x11f
	.long	0x1b66
	.sleb128 480
	.uleb128 0x9
	.long	.LASF343
	.byte	0x12
	.value	0x120
	.long	0x3e0
	.sleb128 484
	.uleb128 0x9
	.long	.LASF344
	.byte	0x12
	.value	0x121
	.long	0x3e0
	.sleb128 488
	.uleb128 0x9
	.long	.LASF345
	.byte	0x12
	.value	0x122
	.long	0x36a
	.sleb128 492
	.uleb128 0x9
	.long	.LASF346
	.byte	0x12
	.value	0x123
	.long	0x3e0
	.sleb128 496
	.uleb128 0x9
	.long	.LASF347
	.byte	0x12
	.value	0x124
	.long	0x1d3b
	.sleb128 504
	.uleb128 0x9
	.long	.LASF329
	.byte	0x12
	.value	0x125
	.long	0xd88
	.sleb128 584
	.uleb128 0x9
	.long	.LASF228
	.byte	0x12
	.value	0x126
	.long	0x3204
	.sleb128 592
	.uleb128 0x9
	.long	.LASF229
	.byte	0x12
	.value	0x126
	.long	0x3204
	.sleb128 744
	.uleb128 0x9
	.long	.LASF348
	.byte	0x12
	.value	0x127
	.long	0x1ff7
	.sleb128 896
	.uleb128 0x9
	.long	.LASF349
	.byte	0x12
	.value	0x128
	.long	0x1ff7
	.sleb128 904
	.uleb128 0x9
	.long	.LASF350
	.byte	0x12
	.value	0x129
	.long	0x3e0
	.sleb128 912
	.uleb128 0x9
	.long	.LASF351
	.byte	0x12
	.value	0x12a
	.long	0x325d
	.sleb128 916
	.uleb128 0x9
	.long	.LASF352
	.byte	0x12
	.value	0x12b
	.long	0x3e0
	.sleb128 956
	.uleb128 0x9
	.long	.LASF353
	.byte	0x12
	.value	0x12b
	.long	0x3e0
	.sleb128 960
	.uleb128 0x9
	.long	.LASF354
	.byte	0x12
	.value	0x12b
	.long	0x3e0
	.sleb128 964
	.uleb128 0x1f
	.string	"spd"
	.byte	0x12
	.value	0x12b
	.long	0x3e0
	.sleb128 968
	.uleb128 0x9
	.long	.LASF355
	.byte	0x12
	.value	0x12c
	.long	0x3e0
	.sleb128 972
	.uleb128 0x9
	.long	.LASF356
	.byte	0x12
	.value	0x12c
	.long	0x3e0
	.sleb128 976
	.uleb128 0x9
	.long	.LASF357
	.byte	0x12
	.value	0x12c
	.long	0x3e0
	.sleb128 980
	.uleb128 0x9
	.long	.LASF358
	.byte	0x12
	.value	0x12c
	.long	0x3e0
	.sleb128 984
	.uleb128 0x9
	.long	.LASF359
	.byte	0x12
	.value	0x12d
	.long	0x3e0
	.sleb128 988
	.uleb128 0x9
	.long	.LASF360
	.byte	0x12
	.value	0x12e
	.long	0x3e0
	.sleb128 992
	.uleb128 0x9
	.long	.LASF361
	.byte	0x12
	.value	0x12e
	.long	0x3e0
	.sleb128 996
	.uleb128 0x9
	.long	.LASF362
	.byte	0x12
	.value	0x130
	.long	0x1913
	.sleb128 1000
	.uleb128 0x9
	.long	.LASF232
	.byte	0x12
	.value	0x132
	.long	0x4b
	.sleb128 1008
	.uleb128 0x9
	.long	.LASF363
	.byte	0x12
	.value	0x133
	.long	0x120
	.sleb128 1016
	.byte	0x0
	.uleb128 0x10
	.byte	0x4
	.byte	0x11
	.byte	0x87
	.long	0x1b66
	.uleb128 0x11
	.long	.LASF364
	.sleb128 0
	.uleb128 0x11
	.long	.LASF365
	.sleb128 1
	.uleb128 0x11
	.long	.LASF366
	.sleb128 2
	.uleb128 0x11
	.long	.LASF367
	.sleb128 3
	.uleb128 0x11
	.long	.LASF368
	.sleb128 4
	.uleb128 0x11
	.long	.LASF369
	.sleb128 5
	.byte	0x0
	.uleb128 0x3
	.long	.LASF370
	.byte	0x11
	.byte	0x87
	.long	0x1b39
	.uleb128 0x10
	.byte	0x4
	.byte	0x11
	.byte	0xa2
	.long	0x1b8c
	.uleb128 0x11
	.long	.LASF371
	.sleb128 0
	.uleb128 0x11
	.long	.LASF372
	.sleb128 1
	.uleb128 0x11
	.long	.LASF373
	.sleb128 2
	.byte	0x0
	.uleb128 0x3
	.long	.LASF374
	.byte	0x11
	.byte	0xa2
	.long	0x1b71
	.uleb128 0x10
	.byte	0x4
	.byte	0x11
	.byte	0xfb
	.long	0x1bb8
	.uleb128 0x11
	.long	.LASF375
	.sleb128 0
	.uleb128 0x11
	.long	.LASF376
	.sleb128 1
	.uleb128 0x11
	.long	.LASF377
	.sleb128 2
	.uleb128 0x11
	.long	.LASF378
	.sleb128 3
	.byte	0x0
	.uleb128 0x3
	.long	.LASF379
	.byte	0x11
	.byte	0xfb
	.long	0x1b97
	.uleb128 0x12
	.byte	0x4
	.byte	0x11
	.value	0x1a1
	.long	0x1bd9
	.uleb128 0x11
	.long	.LASF380
	.sleb128 1
	.uleb128 0x11
	.long	.LASF381
	.sleb128 0
	.byte	0x0
	.uleb128 0x6
	.long	.LASF382
	.byte	0x11
	.value	0x1a1
	.long	0x1bc3
	.uleb128 0x12
	.byte	0x4
	.byte	0x11
	.value	0x1b1
	.long	0x1c79
	.uleb128 0x11
	.long	.LASF383
	.sleb128 0
	.uleb128 0x11
	.long	.LASF384
	.sleb128 1
	.uleb128 0x11
	.long	.LASF385
	.sleb128 2
	.uleb128 0x11
	.long	.LASF386
	.sleb128 3
	.uleb128 0x11
	.long	.LASF387
	.sleb128 4
	.uleb128 0x11
	.long	.LASF388
	.sleb128 5
	.uleb128 0x11
	.long	.LASF389
	.sleb128 6
	.uleb128 0x11
	.long	.LASF390
	.sleb128 7
	.uleb128 0x11
	.long	.LASF391
	.sleb128 8
	.uleb128 0x11
	.long	.LASF392
	.sleb128 9
	.uleb128 0x11
	.long	.LASF393
	.sleb128 10
	.uleb128 0x11
	.long	.LASF394
	.sleb128 11
	.uleb128 0x11
	.long	.LASF395
	.sleb128 12
	.uleb128 0x11
	.long	.LASF396
	.sleb128 13
	.uleb128 0x11
	.long	.LASF397
	.sleb128 14
	.uleb128 0x11
	.long	.LASF398
	.sleb128 15
	.uleb128 0x11
	.long	.LASF399
	.sleb128 16
	.uleb128 0x11
	.long	.LASF400
	.sleb128 17
	.uleb128 0x11
	.long	.LASF401
	.sleb128 18
	.uleb128 0x11
	.long	.LASF402
	.sleb128 19
	.uleb128 0x11
	.long	.LASF403
	.sleb128 20
	.uleb128 0x11
	.long	.LASF404
	.sleb128 21
	.uleb128 0x11
	.long	.LASF405
	.sleb128 22
	.byte	0x0
	.uleb128 0x6
	.long	.LASF406
	.byte	0x11
	.value	0x1bf
	.long	0x1be5
	.uleb128 0x12
	.byte	0x4
	.byte	0x11
	.value	0x1f4
	.long	0x1ca1
	.uleb128 0x11
	.long	.LASF407
	.sleb128 0
	.uleb128 0x11
	.long	.LASF408
	.sleb128 1
	.uleb128 0x11
	.long	.LASF409
	.sleb128 2
	.byte	0x0
	.uleb128 0x6
	.long	.LASF410
	.byte	0x11
	.value	0x1f4
	.long	0x1c85
	.uleb128 0x21
	.byte	0x50
	.byte	0x11
	.value	0x21b
	.long	0x1d3b
	.uleb128 0x9
	.long	.LASF411
	.byte	0x11
	.value	0x21c
	.long	0x3b2
	.sleb128 0
	.uleb128 0x9
	.long	.LASF412
	.byte	0x11
	.value	0x21d
	.long	0x3b2
	.sleb128 8
	.uleb128 0x9
	.long	.LASF413
	.byte	0x11
	.value	0x21d
	.long	0x3b2
	.sleb128 16
	.uleb128 0x9
	.long	.LASF414
	.byte	0x11
	.value	0x21d
	.long	0x3b2
	.sleb128 24
	.uleb128 0x9
	.long	.LASF415
	.byte	0x11
	.value	0x21e
	.long	0x3b2
	.sleb128 32
	.uleb128 0x9
	.long	.LASF416
	.byte	0x11
	.value	0x21f
	.long	0x3b2
	.sleb128 40
	.uleb128 0x9
	.long	.LASF417
	.byte	0x11
	.value	0x220
	.long	0x3b2
	.sleb128 48
	.uleb128 0x9
	.long	.LASF418
	.byte	0x11
	.value	0x221
	.long	0x3b2
	.sleb128 56
	.uleb128 0x9
	.long	.LASF419
	.byte	0x11
	.value	0x221
	.long	0x3b2
	.sleb128 64
	.uleb128 0x9
	.long	.LASF420
	.byte	0x11
	.value	0x222
	.long	0x3b2
	.sleb128 72
	.byte	0x0
	.uleb128 0x6
	.long	.LASF421
	.byte	0x11
	.value	0x223
	.long	0x1cad
	.uleb128 0x12
	.byte	0x4
	.byte	0x11
	.value	0x22f
	.long	0x1d63
	.uleb128 0x11
	.long	.LASF422
	.sleb128 1
	.uleb128 0x11
	.long	.LASF423
	.sleb128 2
	.uleb128 0x11
	.long	.LASF424
	.sleb128 3
	.byte	0x0
	.uleb128 0x6
	.long	.LASF425
	.byte	0x11
	.value	0x22f
	.long	0x1d47
	.uleb128 0x21
	.byte	0x58
	.byte	0x11
	.value	0x48c
	.long	0x1e0a
	.uleb128 0x9
	.long	.LASF426
	.byte	0x11
	.value	0x48d
	.long	0x39c
	.sleb128 0
	.uleb128 0x9
	.long	.LASF427
	.byte	0x11
	.value	0x48e
	.long	0x39c
	.sleb128 8
	.uleb128 0x1f
	.string	"dt"
	.byte	0x11
	.value	0x48f
	.long	0x39c
	.sleb128 16
	.uleb128 0x9
	.long	.LASF428
	.byte	0x11
	.value	0x490
	.long	0x39c
	.sleb128 24
	.uleb128 0x9
	.long	.LASF429
	.byte	0x11
	.value	0x491
	.long	0x39c
	.sleb128 32
	.uleb128 0x9
	.long	.LASF430
	.byte	0x11
	.value	0x492
	.long	0x39c
	.sleb128 40
	.uleb128 0x9
	.long	.LASF431
	.byte	0x11
	.value	0x493
	.long	0x39c
	.sleb128 48
	.uleb128 0x9
	.long	.LASF432
	.byte	0x11
	.value	0x494
	.long	0x39c
	.sleb128 56
	.uleb128 0x9
	.long	.LASF433
	.byte	0x11
	.value	0x496
	.long	0x39c
	.sleb128 64
	.uleb128 0x9
	.long	.LASF434
	.byte	0x11
	.value	0x497
	.long	0x39c
	.sleb128 72
	.uleb128 0x9
	.long	.LASF435
	.byte	0x11
	.value	0x498
	.long	0x39c
	.sleb128 80
	.byte	0x0
	.uleb128 0x6
	.long	.LASF436
	.byte	0x11
	.value	0x499
	.long	0x1d6f
	.uleb128 0x12
	.byte	0x4
	.byte	0x11
	.value	0x4be
	.long	0x1e5e
	.uleb128 0x11
	.long	.LASF437
	.sleb128 1
	.uleb128 0x11
	.long	.LASF438
	.sleb128 2
	.uleb128 0x11
	.long	.LASF439
	.sleb128 3
	.uleb128 0x11
	.long	.LASF440
	.sleb128 4
	.uleb128 0x11
	.long	.LASF441
	.sleb128 8
	.uleb128 0x11
	.long	.LASF442
	.sleb128 12
	.uleb128 0x11
	.long	.LASF443
	.sleb128 16
	.uleb128 0x11
	.long	.LASF444
	.sleb128 32
	.uleb128 0x11
	.long	.LASF445
	.sleb128 64
	.uleb128 0x11
	.long	.LASF446
	.sleb128 128
	.byte	0x0
	.uleb128 0x6
	.long	.LASF447
	.byte	0x11
	.value	0x4c1
	.long	0x1e16
	.uleb128 0x6
	.long	.LASF448
	.byte	0x11
	.value	0x515
	.long	0x1e76
	.uleb128 0x4
	.byte	0x8
	.long	0x1e7c
	.uleb128 0x20
	.long	.LASF449
	.value	0x280
	.byte	0x12
	.value	0x178
	.long	0x1ff7
	.uleb128 0x1f
	.string	"hdr"
	.byte	0x12
	.value	0x179
	.long	0xa00
	.sleb128 0
	.uleb128 0x1f
	.string	"ops"
	.byte	0x12
	.value	0x179
	.long	0x45
	.sleb128 448
	.uleb128 0x1f
	.string	"M"
	.byte	0x12
	.value	0x17a
	.long	0x36a
	.sleb128 456
	.uleb128 0x1f
	.string	"N"
	.byte	0x12
	.value	0x17a
	.long	0x36a
	.sleb128 460
	.uleb128 0x1f
	.string	"m"
	.byte	0x12
	.value	0x17a
	.long	0x36a
	.sleb128 464
	.uleb128 0x9
	.long	.LASF244
	.byte	0x12
	.value	0x17b
	.long	0x36a
	.sleb128 468
	.uleb128 0x9
	.long	.LASF450
	.byte	0x12
	.value	0x17c
	.long	0x36a
	.sleb128 472
	.uleb128 0x9
	.long	.LASF451
	.byte	0x12
	.value	0x17d
	.long	0x991
	.sleb128 480
	.uleb128 0x9
	.long	.LASF452
	.byte	0x12
	.value	0x17e
	.long	0x997
	.sleb128 488
	.uleb128 0x9
	.long	.LASF453
	.byte	0x12
	.value	0x17f
	.long	0x991
	.sleb128 496
	.uleb128 0x9
	.long	.LASF454
	.byte	0x12
	.value	0x180
	.long	0x997
	.sleb128 504
	.uleb128 0x9
	.long	.LASF455
	.byte	0x12
	.value	0x181
	.long	0x997
	.sleb128 512
	.uleb128 0x9
	.long	.LASF456
	.byte	0x12
	.value	0x182
	.long	0x39c
	.sleb128 520
	.uleb128 0x9
	.long	.LASF457
	.byte	0x12
	.value	0x183
	.long	0x39c
	.sleb128 528
	.uleb128 0x1f
	.string	"w1"
	.byte	0x12
	.value	0x184
	.long	0xece
	.sleb128 536
	.uleb128 0x1f
	.string	"w2"
	.byte	0x12
	.value	0x184
	.long	0xece
	.sleb128 544
	.uleb128 0x1f
	.string	"w3"
	.byte	0x12
	.value	0x184
	.long	0xece
	.sleb128 552
	.uleb128 0x1f
	.string	"f"
	.byte	0x12
	.value	0x185
	.long	0x768
	.sleb128 560
	.uleb128 0x9
	.long	.LASF458
	.byte	0x12
	.value	0x186
	.long	0x4b
	.sleb128 568
	.uleb128 0x9
	.long	.LASF459
	.byte	0x12
	.value	0x187
	.long	0x997
	.sleb128 576
	.uleb128 0x9
	.long	.LASF460
	.byte	0x12
	.value	0x188
	.long	0xece
	.sleb128 584
	.uleb128 0x1f
	.string	"F"
	.byte	0x12
	.value	0x189
	.long	0xece
	.sleb128 592
	.uleb128 0x9
	.long	.LASF461
	.byte	0x12
	.value	0x18a
	.long	0x36a
	.sleb128 600
	.uleb128 0x9
	.long	.LASF462
	.byte	0x12
	.value	0x18b
	.long	0x333
	.sleb128 608
	.uleb128 0x9
	.long	.LASF223
	.byte	0x12
	.value	0x18c
	.long	0xe36
	.sleb128 616
	.uleb128 0x9
	.long	.LASF463
	.byte	0x12
	.value	0x18e
	.long	0x4b
	.sleb128 624
	.uleb128 0x9
	.long	.LASF464
	.byte	0x12
	.value	0x18e
	.long	0x4b
	.sleb128 632
	.byte	0x0
	.uleb128 0x6
	.long	.LASF465
	.byte	0x11
	.value	0x65e
	.long	0x2003
	.uleb128 0x4
	.byte	0x8
	.long	0x2009
	.uleb128 0x20
	.long	.LASF466
	.value	0x1f8
	.byte	0x12
	.value	0x194
	.long	0x2093
	.uleb128 0x1f
	.string	"hdr"
	.byte	0x12
	.value	0x195
	.long	0xa00
	.sleb128 0
	.uleb128 0x1f
	.string	"ops"
	.byte	0x12
	.value	0x195
	.long	0x45
	.sleb128 448
	.uleb128 0x9
	.long	.LASF467
	.byte	0x12
	.value	0x196
	.long	0x3e0
	.sleb128 456
	.uleb128 0x1f
	.string	"n"
	.byte	0x12
	.value	0x197
	.long	0x36a
	.sleb128 460
	.uleb128 0x9
	.long	.LASF468
	.byte	0x12
	.value	0x198
	.long	0x13b7
	.sleb128 464
	.uleb128 0x9
	.long	.LASF469
	.byte	0x12
	.value	0x199
	.long	0x9a3
	.sleb128 472
	.uleb128 0x1f
	.string	"vec"
	.byte	0x12
	.value	0x19a
	.long	0xece
	.sleb128 480
	.uleb128 0x9
	.long	.LASF470
	.byte	0x12
	.value	0x19b
	.long	0x32d8
	.sleb128 488
	.uleb128 0x9
	.long	.LASF471
	.byte	0x12
	.value	0x19c
	.long	0x4b
	.sleb128 496
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2099
	.uleb128 0x13
	.long	.LASF472
	.value	0x420
	.byte	0x12
	.byte	0x11
	.long	0x2753
	.uleb128 0xc
	.long	.LASF269
	.byte	0x12
	.byte	0x13
	.long	0x2781
	.sleb128 0
	.uleb128 0xc
	.long	.LASF473
	.byte	0x12
	.byte	0x14
	.long	0x27ab
	.sleb128 8
	.uleb128 0xc
	.long	.LASF474
	.byte	0x12
	.byte	0x15
	.long	0x27ab
	.sleb128 16
	.uleb128 0xc
	.long	.LASF475
	.byte	0x12
	.byte	0x16
	.long	0x27cb
	.sleb128 24
	.uleb128 0xc
	.long	.LASF476
	.byte	0x12
	.byte	0x17
	.long	0x27f0
	.sleb128 32
	.uleb128 0xc
	.long	.LASF477
	.byte	0x12
	.byte	0x19
	.long	0x27cb
	.sleb128 40
	.uleb128 0xc
	.long	.LASF478
	.byte	0x12
	.byte	0x1a
	.long	0x27f0
	.sleb128 48
	.uleb128 0xc
	.long	.LASF479
	.byte	0x12
	.byte	0x1b
	.long	0x27cb
	.sleb128 56
	.uleb128 0xc
	.long	.LASF480
	.byte	0x12
	.byte	0x1c
	.long	0x27f0
	.sleb128 64
	.uleb128 0xc
	.long	.LASF481
	.byte	0x12
	.byte	0x1d
	.long	0x27cb
	.sleb128 72
	.uleb128 0xc
	.long	.LASF482
	.byte	0x12
	.byte	0x1f
	.long	0x27f0
	.sleb128 80
	.uleb128 0xc
	.long	.LASF483
	.byte	0x12
	.byte	0x20
	.long	0x2820
	.sleb128 88
	.uleb128 0xc
	.long	.LASF484
	.byte	0x12
	.byte	0x21
	.long	0x2840
	.sleb128 96
	.uleb128 0x14
	.string	"sor"
	.byte	0x12
	.byte	0x22
	.long	0x2879
	.sleb128 104
	.uleb128 0xc
	.long	.LASF485
	.byte	0x12
	.byte	0x23
	.long	0x289f
	.sleb128 112
	.uleb128 0xc
	.long	.LASF486
	.byte	0x12
	.byte	0x25
	.long	0x28c5
	.sleb128 120
	.uleb128 0xc
	.long	.LASF487
	.byte	0x12
	.byte	0x26
	.long	0x28eb
	.sleb128 128
	.uleb128 0xc
	.long	.LASF488
	.byte	0x12
	.byte	0x27
	.long	0x2906
	.sleb128 136
	.uleb128 0xc
	.long	.LASF489
	.byte	0x12
	.byte	0x28
	.long	0x27cb
	.sleb128 144
	.uleb128 0xc
	.long	.LASF255
	.byte	0x12
	.byte	0x29
	.long	0x2926
	.sleb128 152
	.uleb128 0xc
	.long	.LASF270
	.byte	0x12
	.byte	0x2b
	.long	0x2941
	.sleb128 160
	.uleb128 0xc
	.long	.LASF271
	.byte	0x12
	.byte	0x2c
	.long	0x2941
	.sleb128 168
	.uleb128 0xc
	.long	.LASF277
	.byte	0x12
	.byte	0x2d
	.long	0x2961
	.sleb128 176
	.uleb128 0xc
	.long	.LASF490
	.byte	0x12
	.byte	0x2e
	.long	0x2977
	.sleb128 184
	.uleb128 0xc
	.long	.LASF491
	.byte	0x12
	.byte	0x30
	.long	0x29a6
	.sleb128 192
	.uleb128 0xc
	.long	.LASF492
	.byte	0x12
	.byte	0x31
	.long	0x29d0
	.sleb128 200
	.uleb128 0xc
	.long	.LASF493
	.byte	0x12
	.byte	0x32
	.long	0x29f0
	.sleb128 208
	.uleb128 0xc
	.long	.LASF494
	.byte	0x12
	.byte	0x33
	.long	0x2a15
	.sleb128 216
	.uleb128 0xc
	.long	.LASF495
	.byte	0x12
	.byte	0x34
	.long	0x29f0
	.sleb128 224
	.uleb128 0xc
	.long	.LASF496
	.byte	0x12
	.byte	0x36
	.long	0x2977
	.sleb128 232
	.uleb128 0xc
	.long	.LASF497
	.byte	0x12
	.byte	0x37
	.long	0x29d0
	.sleb128 240
	.uleb128 0xc
	.long	.LASF498
	.byte	0x12
	.byte	0x38
	.long	0x2a15
	.sleb128 248
	.uleb128 0xc
	.long	.LASF272
	.byte	0x12
	.byte	0x39
	.long	0x2a30
	.sleb128 256
	.uleb128 0xc
	.long	.LASF275
	.byte	0x12
	.byte	0x3a
	.long	0x2a30
	.sleb128 264
	.uleb128 0xc
	.long	.LASF251
	.byte	0x12
	.byte	0x3c
	.long	0x2a50
	.sleb128 272
	.uleb128 0xc
	.long	.LASF499
	.byte	0x12
	.byte	0x3d
	.long	0x27cb
	.sleb128 280
	.uleb128 0xc
	.long	.LASF500
	.byte	0x12
	.byte	0x3e
	.long	0x27cb
	.sleb128 288
	.uleb128 0xc
	.long	.LASF501
	.byte	0x12
	.byte	0x3f
	.long	0x2820
	.sleb128 296
	.uleb128 0xc
	.long	.LASF502
	.byte	0x12
	.byte	0x40
	.long	0x2840
	.sleb128 304
	.uleb128 0xc
	.long	.LASF261
	.byte	0x12
	.byte	0x42
	.long	0x2a75
	.sleb128 312
	.uleb128 0xc
	.long	.LASF503
	.byte	0x12
	.byte	0x43
	.long	0x2ab5
	.sleb128 320
	.uleb128 0xc
	.long	.LASF504
	.byte	0x12
	.byte	0x44
	.long	0x2ada
	.sleb128 328
	.uleb128 0xc
	.long	.LASF297
	.byte	0x12
	.byte	0x45
	.long	0x2b09
	.sleb128 336
	.uleb128 0xc
	.long	.LASF259
	.byte	0x12
	.byte	0x46
	.long	0x2b29
	.sleb128 344
	.uleb128 0xc
	.long	.LASF505
	.byte	0x12
	.byte	0x48
	.long	0x2b49
	.sleb128 352
	.uleb128 0xc
	.long	.LASF258
	.byte	0x12
	.byte	0x49
	.long	0x2b64
	.sleb128 360
	.uleb128 0xc
	.long	.LASF299
	.byte	0x12
	.byte	0x4a
	.long	0x2b64
	.sleb128 368
	.uleb128 0xc
	.long	.LASF506
	.byte	0x12
	.byte	0x4b
	.long	0x2b84
	.sleb128 376
	.uleb128 0xc
	.long	.LASF507
	.byte	0x12
	.byte	0x4c
	.long	0x29a6
	.sleb128 384
	.uleb128 0xc
	.long	.LASF508
	.byte	0x12
	.byte	0x4e
	.long	0x2b9f
	.sleb128 392
	.uleb128 0xc
	.long	.LASF509
	.byte	0x12
	.byte	0x4f
	.long	0x2bd8
	.sleb128 400
	.uleb128 0xc
	.long	.LASF510
	.byte	0x12
	.byte	0x50
	.long	0x2bd8
	.sleb128 408
	.uleb128 0xc
	.long	.LASF511
	.byte	0x12
	.byte	0x51
	.long	0x2bd8
	.sleb128 416
	.uleb128 0xc
	.long	.LASF512
	.byte	0x12
	.byte	0x52
	.long	0x2bd8
	.sleb128 424
	.uleb128 0xc
	.long	.LASF513
	.byte	0x12
	.byte	0x54
	.long	0x2bf8
	.sleb128 432
	.uleb128 0xc
	.long	.LASF514
	.byte	0x12
	.byte	0x55
	.long	0x2c28
	.sleb128 440
	.uleb128 0xc
	.long	.LASF515
	.byte	0x12
	.byte	0x56
	.long	0x2977
	.sleb128 448
	.uleb128 0xc
	.long	.LASF516
	.byte	0x12
	.byte	0x57
	.long	0x2c4d
	.sleb128 456
	.uleb128 0xc
	.long	.LASF278
	.byte	0x12
	.byte	0x58
	.long	0x2781
	.sleb128 464
	.uleb128 0xc
	.long	.LASF517
	.byte	0x12
	.byte	0x5a
	.long	0x2c77
	.sleb128 472
	.uleb128 0xc
	.long	.LASF152
	.byte	0x12
	.byte	0x5b
	.long	0x2977
	.sleb128 480
	.uleb128 0xc
	.long	.LASF151
	.byte	0x12
	.byte	0x5c
	.long	0x2c92
	.sleb128 488
	.uleb128 0xc
	.long	.LASF518
	.byte	0x12
	.byte	0x5d
	.long	0x2cb7
	.sleb128 496
	.uleb128 0xc
	.long	.LASF519
	.byte	0x12
	.byte	0x5e
	.long	0x2cd2
	.sleb128 504
	.uleb128 0xc
	.long	.LASF520
	.byte	0x12
	.byte	0x60
	.long	0x27cb
	.sleb128 512
	.uleb128 0xc
	.long	.LASF521
	.byte	0x12
	.byte	0x61
	.long	0x27cb
	.sleb128 520
	.uleb128 0xc
	.long	.LASF289
	.byte	0x12
	.byte	0x62
	.long	0x2cf2
	.sleb128 528
	.uleb128 0xc
	.long	.LASF290
	.byte	0x12
	.byte	0x63
	.long	0x2781
	.sleb128 536
	.uleb128 0xc
	.long	.LASF522
	.byte	0x12
	.byte	0x64
	.long	0x29a6
	.sleb128 544
	.uleb128 0xc
	.long	.LASF523
	.byte	0x12
	.byte	0x66
	.long	0x2b49
	.sleb128 552
	.uleb128 0xc
	.long	.LASF524
	.byte	0x12
	.byte	0x67
	.long	0x2b49
	.sleb128 560
	.uleb128 0xc
	.long	.LASF525
	.byte	0x12
	.byte	0x68
	.long	0x2cb7
	.sleb128 568
	.uleb128 0xc
	.long	.LASF526
	.byte	0x12
	.byte	0x69
	.long	0x2d0d
	.sleb128 576
	.uleb128 0xc
	.long	.LASF527
	.byte	0x12
	.byte	0x6a
	.long	0x2d28
	.sleb128 584
	.uleb128 0xc
	.long	.LASF528
	.byte	0x12
	.byte	0x6c
	.long	0x2d48
	.sleb128 592
	.uleb128 0xc
	.long	.LASF529
	.byte	0x12
	.byte	0x6d
	.long	0x2d78
	.sleb128 600
	.uleb128 0xc
	.long	.LASF292
	.byte	0x12
	.byte	0x6e
	.long	0x2977
	.sleb128 608
	.uleb128 0xc
	.long	.LASF530
	.byte	0x12
	.byte	0x6f
	.long	0x27cb
	.sleb128 616
	.uleb128 0xc
	.long	.LASF531
	.byte	0x12
	.byte	0x70
	.long	0x27cb
	.sleb128 624
	.uleb128 0xc
	.long	.LASF532
	.byte	0x12
	.byte	0x72
	.long	0x2d93
	.sleb128 632
	.uleb128 0xc
	.long	.LASF533
	.byte	0x12
	.byte	0x73
	.long	0x2db3
	.sleb128 640
	.uleb128 0xc
	.long	.LASF534
	.byte	0x12
	.byte	0x74
	.long	0x2db3
	.sleb128 648
	.uleb128 0xc
	.long	.LASF535
	.byte	0x12
	.byte	0x75
	.long	0x2dd8
	.sleb128 656
	.uleb128 0xc
	.long	.LASF286
	.byte	0x12
	.byte	0x76
	.long	0x2c92
	.sleb128 664
	.uleb128 0xc
	.long	.LASF536
	.byte	0x12
	.byte	0x78
	.long	0x2df8
	.sleb128 672
	.uleb128 0xc
	.long	.LASF537
	.byte	0x12
	.byte	0x79
	.long	0x2df8
	.sleb128 680
	.uleb128 0xc
	.long	.LASF538
	.byte	0x12
	.byte	0x7a
	.long	0x2e13
	.sleb128 688
	.uleb128 0xc
	.long	.LASF539
	.byte	0x12
	.byte	0x7b
	.long	0x2781
	.sleb128 696
	.uleb128 0xc
	.long	.LASF540
	.byte	0x12
	.byte	0x7c
	.long	0x2e33
	.sleb128 704
	.uleb128 0xc
	.long	.LASF541
	.byte	0x12
	.byte	0x7e
	.long	0x2e5d
	.sleb128 712
	.uleb128 0xc
	.long	.LASF542
	.byte	0x12
	.byte	0x7f
	.long	0x2e82
	.sleb128 720
	.uleb128 0xc
	.long	.LASF543
	.byte	0x12
	.byte	0x80
	.long	0x2ea2
	.sleb128 728
	.uleb128 0xc
	.long	.LASF544
	.byte	0x12
	.byte	0x81
	.long	0x2e5d
	.sleb128 736
	.uleb128 0xc
	.long	.LASF545
	.byte	0x12
	.byte	0x82
	.long	0x2e82
	.sleb128 744
	.uleb128 0xc
	.long	.LASF546
	.byte	0x12
	.byte	0x84
	.long	0x2ea2
	.sleb128 752
	.uleb128 0xc
	.long	.LASF547
	.byte	0x12
	.byte	0x85
	.long	0x2e5d
	.sleb128 760
	.uleb128 0xc
	.long	.LASF548
	.byte	0x12
	.byte	0x86
	.long	0x2e82
	.sleb128 768
	.uleb128 0xc
	.long	.LASF549
	.byte	0x12
	.byte	0x87
	.long	0x2ea2
	.sleb128 776
	.uleb128 0xc
	.long	.LASF550
	.byte	0x12
	.byte	0x88
	.long	0x2e82
	.sleb128 784
	.uleb128 0xc
	.long	.LASF551
	.byte	0x12
	.byte	0x8a
	.long	0x2ea2
	.sleb128 792
	.uleb128 0xc
	.long	.LASF552
	.byte	0x12
	.byte	0x8b
	.long	0x2e82
	.sleb128 800
	.uleb128 0xc
	.long	.LASF553
	.byte	0x12
	.byte	0x8c
	.long	0x2ea2
	.sleb128 808
	.uleb128 0xc
	.long	.LASF288
	.byte	0x12
	.byte	0x8d
	.long	0x2977
	.sleb128 816
	.uleb128 0xc
	.long	.LASF554
	.byte	0x12
	.byte	0x8e
	.long	0x2ecc
	.sleb128 824
	.uleb128 0xc
	.long	.LASF555
	.byte	0x12
	.byte	0x90
	.long	0x2eec
	.sleb128 832
	.uleb128 0xc
	.long	.LASF556
	.byte	0x12
	.byte	0x91
	.long	0x2977
	.sleb128 840
	.uleb128 0xc
	.long	.LASF557
	.byte	0x12
	.byte	0x92
	.long	0x2977
	.sleb128 848
	.uleb128 0xc
	.long	.LASF558
	.byte	0x12
	.byte	0x93
	.long	0x2977
	.sleb128 856
	.uleb128 0xc
	.long	.LASF559
	.byte	0x12
	.byte	0x94
	.long	0x2977
	.sleb128 864
	.uleb128 0xc
	.long	.LASF560
	.byte	0x12
	.byte	0x96
	.long	0x2ea2
	.sleb128 872
	.uleb128 0xc
	.long	.LASF561
	.byte	0x12
	.byte	0x97
	.long	0x2f1b
	.sleb128 880
	.uleb128 0xc
	.long	.LASF562
	.byte	0x12
	.byte	0x98
	.long	0x2b49
	.sleb128 888
	.uleb128 0xc
	.long	.LASF563
	.byte	0x12
	.byte	0x99
	.long	0x2f3b
	.sleb128 896
	.uleb128 0xc
	.long	.LASF564
	.byte	0x12
	.byte	0x9a
	.long	0x2f5b
	.sleb128 904
	.uleb128 0xc
	.long	.LASF565
	.byte	0x12
	.byte	0x9c
	.long	0x2f76
	.sleb128 912
	.uleb128 0xc
	.long	.LASF300
	.byte	0x12
	.byte	0x9d
	.long	0x2977
	.sleb128 920
	.uleb128 0xc
	.long	.LASF566
	.byte	0x12
	.byte	0x9e
	.long	0x2f9c
	.sleb128 928
	.uleb128 0xc
	.long	.LASF567
	.byte	0x12
	.byte	0x9f
	.long	0x2c4d
	.sleb128 936
	.uleb128 0xc
	.long	.LASF568
	.byte	0x12
	.byte	0xa0
	.long	0x2c4d
	.sleb128 944
	.uleb128 0xc
	.long	.LASF569
	.byte	0x12
	.byte	0xa2
	.long	0x27cb
	.sleb128 952
	.uleb128 0xc
	.long	.LASF570
	.byte	0x12
	.byte	0xa3
	.long	0x289f
	.sleb128 960
	.uleb128 0xc
	.long	.LASF571
	.byte	0x12
	.byte	0xa4
	.long	0x27cb
	.sleb128 968
	.uleb128 0xc
	.long	.LASF572
	.byte	0x12
	.byte	0xa5
	.long	0x27f0
	.sleb128 976
	.uleb128 0xc
	.long	.LASF573
	.byte	0x12
	.byte	0xa6
	.long	0x2fbc
	.sleb128 984
	.uleb128 0xc
	.long	.LASF574
	.byte	0x12
	.byte	0xa8
	.long	0x2d93
	.sleb128 992
	.uleb128 0xc
	.long	.LASF575
	.byte	0x12
	.byte	0xa9
	.long	0x2926
	.sleb128 1000
	.uleb128 0xc
	.long	.LASF576
	.byte	0x12
	.byte	0xaa
	.long	0x2a30
	.sleb128 1008
	.uleb128 0xc
	.long	.LASF577
	.byte	0x12
	.byte	0xab
	.long	0x27f0
	.sleb128 1016
	.uleb128 0xc
	.long	.LASF578
	.byte	0x12
	.byte	0xac
	.long	0x2ab5
	.sleb128 1024
	.uleb128 0xc
	.long	.LASF579
	.byte	0x12
	.byte	0xae
	.long	0x2fe6
	.sleb128 1032
	.uleb128 0xc
	.long	.LASF580
	.byte	0x12
	.byte	0xb0
	.long	0x3010
	.sleb128 1040
	.uleb128 0xc
	.long	.LASF581
	.byte	0x12
	.byte	0xb2
	.long	0x3035
	.sleb128 1048
	.byte	0x0
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2781
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x150e
	.uleb128 0x1d
	.long	0xd88
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2753
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x27ab
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x991
	.uleb128 0x1d
	.long	0x997
	.uleb128 0x1d
	.long	0x9a9
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2787
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x27cb
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x27b1
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x27f0
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x27d1
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2815
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0x2815
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x281b
	.uleb128 0xf
	.long	0x1e0a
	.uleb128 0x4
	.byte	0x8
	.long	0x27f6
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2840
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0x2815
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2826
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2879
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x39c
	.uleb128 0x1d
	.long	0x1e5e
	.uleb128 0x1d
	.long	0x39c
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2846
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2899
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1b8c
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x195a
	.uleb128 0x4
	.byte	0x8
	.long	0x287f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x28bf
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1d63
	.uleb128 0x1d
	.long	0x28bf
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1d3b
	.uleb128 0x4
	.byte	0x8
	.long	0x28a5
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x28e5
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x28e5
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x3e0
	.uleb128 0x4
	.byte	0x8
	.long	0x28cb
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2906
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x28f1
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2926
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xf98
	.uleb128 0x1d
	.long	0x824
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x290c
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2941
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1bd9
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x292c
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2961
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1c79
	.uleb128 0x1d
	.long	0x3e0
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2947
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2977
	.uleb128 0x1d
	.long	0x195a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2967
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x29a6
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xece
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x297d
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x29d0
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0x2815
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x29ac
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x29f0
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x2815
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x29d6
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2a15
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0x2815
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x29f6
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2a30
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x9a9
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2a1b
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2a50
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1ca1
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2a36
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2a75
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x3a7
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1bb8
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2a56
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2aa4
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x2aa4
	.uleb128 0x1d
	.long	0x2aa4
	.uleb128 0x1d
	.long	0x1b8c
	.uleb128 0x1d
	.long	0x2aaf
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2aaa
	.uleb128 0xf
	.long	0xd94
	.uleb128 0x4
	.byte	0x8
	.long	0x2899
	.uleb128 0x4
	.byte	0x8
	.long	0x2a7b
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2ada
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0xeb0
	.uleb128 0x1d
	.long	0x36a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2abb
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2b09
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x9a3
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2ae0
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2b29
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1bb8
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2b0f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2b49
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x991
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2b2f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2b64
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x3a7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2b4f
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2b84
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0xd88
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2b6a
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2b9f
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2b8a
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2bd8
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x3e0
	.uleb128 0x1d
	.long	0x3e0
	.uleb128 0x1d
	.long	0x991
	.uleb128 0x1d
	.long	0x997
	.uleb128 0x1d
	.long	0x997
	.uleb128 0x1d
	.long	0x28e5
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2ba5
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2bf8
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xebc
	.uleb128 0x1d
	.long	0x1e6a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2bde
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2c22
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0xeb6
	.uleb128 0x1d
	.long	0x2c22
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0xebc
	.uleb128 0x4
	.byte	0x8
	.long	0x2bfe
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2c4d
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2c2e
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2c77
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0xd94
	.uleb128 0x1d
	.long	0x1b8c
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2c53
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2c92
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x6cc
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2c7d
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2cb7
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x1b8c
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2c98
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2cd2
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x3e0
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2cbd
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2cf2
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xe10
	.uleb128 0x1d
	.long	0xe10
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2cd8
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2d0d
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xebc
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2cf8
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2d28
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x4b
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2d13
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2d48
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x4b
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2d2e
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2d72
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1e6a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x2d72
	.uleb128 0x1d
	.long	0x4b
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x1bb8
	.uleb128 0x4
	.byte	0x8
	.long	0x2d4e
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2d93
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xeb0
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2d7e
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2db3
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1948
	.uleb128 0x1d
	.long	0x1948
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2d99
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2dd8
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x991
	.uleb128 0x1d
	.long	0x991
	.uleb128 0x1d
	.long	0x991
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2db9
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2df8
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x39c
	.uleb128 0x1d
	.long	0x28e5
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2dde
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2e13
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x28e5
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2dfe
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2e33
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x13b7
	.uleb128 0x1d
	.long	0x13b7
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2e19
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2e5d
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x1b8c
	.uleb128 0x1d
	.long	0x39c
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2e39
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2e82
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x39c
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2e63
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2ea2
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x195a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2e88
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2ecc
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2ea8
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2eec
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x150e
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2ed2
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2f1b
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x34
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x1b8c
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2ef2
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2f3b
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x36a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2f21
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2f5b
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x28e5
	.uleb128 0x1d
	.long	0x991
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2f41
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2f76
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2f61
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2f96
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x991
	.uleb128 0x1d
	.long	0x2f96
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x15b7
	.uleb128 0x4
	.byte	0x8
	.long	0x2f7c
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2fbc
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x34
	.uleb128 0x1d
	.long	0x2899
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2fa2
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x2fe6
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x991
	.uleb128 0x1d
	.long	0x150e
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2fc2
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x3010
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x15b7
	.uleb128 0x1d
	.long	0x36a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x2fec
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x3035
	.uleb128 0x1d
	.long	0x195a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.uleb128 0x1d
	.long	0x36a
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x3016
	.uleb128 0x3
	.long	.LASF582
	.byte	0x12
	.byte	0xd5
	.long	0x3046
	.uleb128 0x4
	.byte	0x8
	.long	0x304c
	.uleb128 0xb
	.long	.LASF583
	.byte	0x38
	.byte	0x12
	.byte	0xd7
	.long	0x30b9
	.uleb128 0xc
	.long	.LASF584
	.byte	0x12
	.byte	0xd8
	.long	0x303b
	.sleb128 0
	.uleb128 0xc
	.long	.LASF585
	.byte	0x12
	.byte	0xd9
	.long	0x9a3
	.sleb128 8
	.uleb128 0x14
	.string	"val"
	.byte	0x12
	.byte	0xd9
	.long	0x9a3
	.sleb128 16
	.uleb128 0x14
	.string	"idx"
	.byte	0x12
	.byte	0xda
	.long	0x991
	.sleb128 24
	.uleb128 0x14
	.string	"idy"
	.byte	0x12
	.byte	0xda
	.long	0x991
	.sleb128 32
	.uleb128 0xc
	.long	.LASF586
	.byte	0x12
	.byte	0xdb
	.long	0x36a
	.sleb128 40
	.uleb128 0xc
	.long	.LASF587
	.byte	0x12
	.byte	0xdc
	.long	0x36a
	.sleb128 44
	.uleb128 0xc
	.long	.LASF588
	.byte	0x12
	.byte	0xdd
	.long	0x36a
	.sleb128 48
	.byte	0x0
	.uleb128 0x1b
	.byte	0x98
	.byte	0x12
	.byte	0xe4
	.long	0x3204
	.uleb128 0xc
	.long	.LASF306
	.byte	0x12
	.byte	0xe5
	.long	0x36a
	.sleb128 0
	.uleb128 0xc
	.long	.LASF307
	.byte	0x12
	.byte	0xe6
	.long	0x36a
	.sleb128 4
	.uleb128 0xc
	.long	.LASF308
	.byte	0x12
	.byte	0xe7
	.long	0x36a
	.sleb128 8
	.uleb128 0x14
	.string	"n"
	.byte	0x12
	.byte	0xe8
	.long	0x36a
	.sleb128 12
	.uleb128 0x14
	.string	"bs"
	.byte	0x12
	.byte	0xe9
	.long	0x36a
	.sleb128 16
	.uleb128 0xc
	.long	.LASF309
	.byte	0x12
	.byte	0xea
	.long	0x36a
	.sleb128 20
	.uleb128 0xc
	.long	.LASF585
	.byte	0x12
	.byte	0xeb
	.long	0x303b
	.sleb128 24
	.uleb128 0xc
	.long	.LASF589
	.byte	0x12
	.byte	0xeb
	.long	0x303b
	.sleb128 32
	.uleb128 0xc
	.long	.LASF78
	.byte	0x12
	.byte	0xed
	.long	0x34
	.sleb128 40
	.uleb128 0xc
	.long	.LASF311
	.byte	0x12
	.byte	0xee
	.long	0x35f
	.sleb128 44
	.uleb128 0xc
	.long	.LASF312
	.byte	0x12
	.byte	0xee
	.long	0x35f
	.sleb128 48
	.uleb128 0xc
	.long	.LASF313
	.byte	0x12
	.byte	0xef
	.long	0x35f
	.sleb128 52
	.uleb128 0xc
	.long	.LASF314
	.byte	0x12
	.byte	0xef
	.long	0x35f
	.sleb128 56
	.uleb128 0xc
	.long	.LASF315
	.byte	0x12
	.byte	0xf0
	.long	0x18e1
	.sleb128 64
	.uleb128 0xc
	.long	.LASF316
	.byte	0x12
	.byte	0xf1
	.long	0x18e1
	.sleb128 72
	.uleb128 0xc
	.long	.LASF317
	.byte	0x12
	.byte	0xf2
	.long	0xc2
	.sleb128 80
	.uleb128 0xc
	.long	.LASF318
	.byte	0x12
	.byte	0xf3
	.long	0x36a
	.sleb128 88
	.uleb128 0xc
	.long	.LASF319
	.byte	0x12
	.byte	0xf3
	.long	0x36a
	.sleb128 92
	.uleb128 0xc
	.long	.LASF320
	.byte	0x12
	.byte	0xf4
	.long	0x9a3
	.sleb128 96
	.uleb128 0xc
	.long	.LASF322
	.byte	0x12
	.byte	0xf5
	.long	0x991
	.sleb128 104
	.uleb128 0xc
	.long	.LASF321
	.byte	0x12
	.byte	0xf6
	.long	0x9a9
	.sleb128 112
	.uleb128 0xc
	.long	.LASF323
	.byte	0x12
	.byte	0xf7
	.long	0x997
	.sleb128 120
	.uleb128 0xc
	.long	.LASF326
	.byte	0x12
	.byte	0xf8
	.long	0x36a
	.sleb128 128
	.uleb128 0xc
	.long	.LASF590
	.byte	0x12
	.byte	0xf9
	.long	0x191e
	.sleb128 136
	.uleb128 0xc
	.long	.LASF591
	.byte	0x12
	.byte	0xfa
	.long	0x3e0
	.sleb128 144
	.uleb128 0xc
	.long	.LASF592
	.byte	0x12
	.byte	0xfb
	.long	0x36a
	.sleb128 148
	.byte	0x0
	.uleb128 0x3
	.long	.LASF593
	.byte	0x12
	.byte	0xfc
	.long	0x30b9
	.uleb128 0x21
	.byte	0x28
	.byte	0x12
	.value	0x10a
	.long	0x324d
	.uleb128 0x1f
	.string	"dim"
	.byte	0x12
	.value	0x10b
	.long	0x36a
	.sleb128 0
	.uleb128 0x9
	.long	.LASF594
	.byte	0x12
	.value	0x10c
	.long	0x324d
	.sleb128 4
	.uleb128 0x9
	.long	.LASF595
	.byte	0x12
	.value	0x10d
	.long	0x324d
	.sleb128 20
	.uleb128 0x1f
	.string	"noc"
	.byte	0x12
	.value	0x10e
	.long	0x3e0
	.sleb128 36
	.byte	0x0
	.uleb128 0xd
	.long	0x36a
	.long	0x325d
	.uleb128 0xe
	.long	0xe0
	.byte	0x3
	.byte	0x0
	.uleb128 0x6
	.long	.LASF596
	.byte	0x12
	.value	0x10f
	.long	0x320f
	.uleb128 0x21
	.byte	0x20
	.byte	0x12
	.value	0x112
	.long	0x32b2
	.uleb128 0x9
	.long	.LASF597
	.byte	0x12
	.value	0x113
	.long	0x3e0
	.sleb128 0
	.uleb128 0x1f
	.string	"use"
	.byte	0x12
	.value	0x114
	.long	0x3e0
	.sleb128 4
	.uleb128 0x9
	.long	.LASF453
	.byte	0x12
	.value	0x115
	.long	0x36a
	.sleb128 8
	.uleb128 0x1f
	.string	"i"
	.byte	0x12
	.value	0x116
	.long	0x991
	.sleb128 16
	.uleb128 0x9
	.long	.LASF598
	.byte	0x12
	.value	0x117
	.long	0x991
	.sleb128 24
	.byte	0x0
	.uleb128 0x6
	.long	.LASF599
	.byte	0x12
	.value	0x118
	.long	0x3269
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x32d8
	.uleb128 0x1d
	.long	0x1ff7
	.uleb128 0x1d
	.long	0xece
	.uleb128 0x1d
	.long	0x4b
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x32be
	.uleb128 0xd
	.long	0x126
	.long	0x32ee
	.uleb128 0xe
	.long	0xe0
	.byte	0xf
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x3be
	.uleb128 0x22
	.value	0x120
	.byte	0x13
	.byte	0x1c
	.long	0x34e6
	.uleb128 0xc
	.long	.LASF600
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 0
	.uleb128 0xc
	.long	.LASF601
	.byte	0x13
	.byte	0x1d
	.long	0x36a
	.sleb128 4
	.uleb128 0xc
	.long	.LASF602
	.byte	0x13
	.byte	0x1d
	.long	0x36a
	.sleb128 8
	.uleb128 0xc
	.long	.LASF603
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 12
	.uleb128 0xc
	.long	.LASF604
	.byte	0x13
	.byte	0x1d
	.long	0x36a
	.sleb128 16
	.uleb128 0xc
	.long	.LASF605
	.byte	0x13
	.byte	0x1d
	.long	0x991
	.sleb128 24
	.uleb128 0xc
	.long	.LASF606
	.byte	0x13
	.byte	0x1d
	.long	0x991
	.sleb128 32
	.uleb128 0xc
	.long	.LASF607
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 40
	.uleb128 0xc
	.long	.LASF309
	.byte	0x13
	.byte	0x1d
	.long	0x36a
	.sleb128 44
	.uleb128 0xc
	.long	.LASF324
	.byte	0x13
	.byte	0x1d
	.long	0x36a
	.sleb128 48
	.uleb128 0xc
	.long	.LASF608
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 52
	.uleb128 0xc
	.long	.LASF609
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 56
	.uleb128 0xc
	.long	.LASF610
	.byte	0x13
	.byte	0x1d
	.long	0x991
	.sleb128 64
	.uleb128 0xc
	.long	.LASF611
	.byte	0x13
	.byte	0x1d
	.long	0x991
	.sleb128 72
	.uleb128 0xc
	.long	.LASF612
	.byte	0x13
	.byte	0x1d
	.long	0x195a
	.sleb128 80
	.uleb128 0xc
	.long	.LASF613
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 88
	.uleb128 0xc
	.long	.LASF614
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 92
	.uleb128 0xc
	.long	.LASF615
	.byte	0x13
	.byte	0x1d
	.long	0x32b2
	.sleb128 96
	.uleb128 0x14
	.string	"nz"
	.byte	0x13
	.byte	0x1d
	.long	0x36a
	.sleb128 128
	.uleb128 0x14
	.string	"i"
	.byte	0x13
	.byte	0x1d
	.long	0x991
	.sleb128 136
	.uleb128 0x14
	.string	"j"
	.byte	0x13
	.byte	0x1d
	.long	0x991
	.sleb128 144
	.uleb128 0xc
	.long	.LASF616
	.byte	0x13
	.byte	0x1d
	.long	0x991
	.sleb128 152
	.uleb128 0xc
	.long	.LASF617
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 160
	.uleb128 0x14
	.string	"a"
	.byte	0x13
	.byte	0x1d
	.long	0x32ee
	.sleb128 168
	.uleb128 0xc
	.long	.LASF618
	.byte	0x13
	.byte	0x1d
	.long	0x9a3
	.sleb128 176
	.uleb128 0x14
	.string	"row"
	.byte	0x13
	.byte	0x1d
	.long	0xd94
	.sleb128 184
	.uleb128 0x14
	.string	"col"
	.byte	0x13
	.byte	0x1d
	.long	0xd94
	.sleb128 192
	.uleb128 0xc
	.long	.LASF619
	.byte	0x13
	.byte	0x1d
	.long	0xd94
	.sleb128 200
	.uleb128 0xc
	.long	.LASF432
	.byte	0x13
	.byte	0x1d
	.long	0x3e0
	.sleb128 208
	.uleb128 0xc
	.long	.LASF89
	.byte	0x13
	.byte	0x1d
	.long	0x195a
	.sleb128 216
	.uleb128 0x14
	.string	"bs2"
	.byte	0x13
	.byte	0x1e
	.long	0x36a
	.sleb128 224
	.uleb128 0x14
	.string	"mbs"
	.byte	0x13
	.byte	0x1e
	.long	0x36a
	.sleb128 228
	.uleb128 0x14
	.string	"nbs"
	.byte	0x13
	.byte	0x1e
	.long	0x36a
	.sleb128 232
	.uleb128 0xc
	.long	.LASF620
	.byte	0x13
	.byte	0x1e
	.long	0x9a3
	.sleb128 240
	.uleb128 0xc
	.long	.LASF621
	.byte	0x13
	.byte	0x1e
	.long	0x9a3
	.sleb128 248
	.uleb128 0xc
	.long	.LASF622
	.byte	0x13
	.byte	0x1e
	.long	0x32ee
	.sleb128 256
	.uleb128 0xc
	.long	.LASF623
	.byte	0x13
	.byte	0x1e
	.long	0x195a
	.sleb128 264
	.uleb128 0xc
	.long	.LASF624
	.byte	0x13
	.byte	0x1e
	.long	0x32ee
	.sleb128 272
	.uleb128 0xc
	.long	.LASF625
	.byte	0x13
	.byte	0x1e
	.long	0x3e0
	.sleb128 280
	.byte	0x0
	.uleb128 0x3
	.long	.LASF626
	.byte	0x13
	.byte	0x1f
	.long	0x32f4
	.uleb128 0x3
	.long	.LASF627
	.byte	0x14
	.byte	0x1e
	.long	0x120
	.uleb128 0x23
	.long	.LASF628
	.byte	0x1
	.byte	0x86
	.byte	0x1
	.long	0x39c
	.quad	.LFB0
	.quad	.LFE0
	.byte	0x1
	.byte	0x9c
	.long	0x352b
	.uleb128 0x24
	.string	"a"
	.byte	0x1
	.byte	0x86
	.long	0x3a7
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.uleb128 0x25
	.long	.LASF629
	.byte	0x2
	.value	0x6e2
	.byte	0x1
	.long	0x33e
	.quad	.LFB6
	.quad	.LFE6
	.byte	0x1
	.byte	0x9c
	.long	0x358a
	.uleb128 0x26
	.string	"a"
	.byte	0x2
	.value	0x6e2
	.long	0x4b
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.uleb128 0x26
	.string	"b"
	.byte	0x2
	.value	0x6e2
	.long	0x710
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x26
	.string	"n"
	.byte	0x2
	.value	0x6e2
	.long	0xd5
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x27
	.long	.LASF632
	.long	0x359a
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.12688
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x359a
	.uleb128 0xe
	.long	0xe0
	.byte	0xb
	.byte	0x0
	.uleb128 0xf
	.long	0x358a
	.uleb128 0x25
	.long	.LASF630
	.byte	0x2
	.value	0x722
	.byte	0x1
	.long	0x33e
	.quad	.LFB7
	.quad	.LFE7
	.byte	0x1
	.byte	0x9c
	.long	0x35dd
	.uleb128 0x26
	.string	"a"
	.byte	0x2
	.value	0x722
	.long	0x4b
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.uleb128 0x26
	.string	"n"
	.byte	0x2
	.value	0x722
	.long	0xd5
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.byte	0x0
	.uleb128 0x25
	.long	.LASF631
	.byte	0x3
	.value	0x115
	.byte	0x1
	.long	0x33e
	.quad	.LFB15
	.quad	.LFE15
	.byte	0x1
	.byte	0x9c
	.long	0x363e
	.uleb128 0x26
	.string	"x"
	.byte	0x3
	.value	0x115
	.long	0xece
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x26
	.string	"a"
	.byte	0x3
	.value	0x115
	.long	0x363e
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF634
	.byte	0x3
	.value	0x117
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x27
	.long	.LASF632
	.long	0x3644
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.15175
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x150e
	.uleb128 0xf
	.long	0x32de
	.uleb128 0x25
	.long	.LASF633
	.byte	0x3
	.value	0x129
	.byte	0x1
	.long	0x33e
	.quad	.LFB16
	.quad	.LFE16
	.byte	0x1
	.byte	0x9c
	.long	0x36aa
	.uleb128 0x26
	.string	"x"
	.byte	0x3
	.value	0x129
	.long	0xece
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x26
	.string	"a"
	.byte	0x3
	.value	0x129
	.long	0x363e
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF634
	.byte	0x3
	.value	0x12b
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x27
	.long	.LASF632
	.long	0x36aa
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.15257
	.byte	0x0
	.uleb128 0xf
	.long	0x323
	.uleb128 0x25
	.long	.LASF635
	.byte	0x3
	.value	0x13c
	.byte	0x1
	.long	0x33e
	.quad	.LFB17
	.quad	.LFE17
	.byte	0x1
	.byte	0x9c
	.long	0x3710
	.uleb128 0x26
	.string	"x"
	.byte	0x3
	.value	0x13c
	.long	0xece
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x26
	.string	"a"
	.byte	0x3
	.value	0x13c
	.long	0x9a9
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF634
	.byte	0x3
	.value	0x13e
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x27
	.long	.LASF632
	.long	0x3710
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.15326
	.byte	0x0
	.uleb128 0xf
	.long	0x358a
	.uleb128 0x25
	.long	.LASF636
	.byte	0x3
	.value	0x150
	.byte	0x1
	.long	0x33e
	.quad	.LFB18
	.quad	.LFE18
	.byte	0x1
	.byte	0x9c
	.long	0x3776
	.uleb128 0x26
	.string	"x"
	.byte	0x3
	.value	0x150
	.long	0xece
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x26
	.string	"a"
	.byte	0x3
	.value	0x150
	.long	0x9a9
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF634
	.byte	0x3
	.value	0x152
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x27
	.long	.LASF632
	.long	0x3776
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.15408
	.byte	0x0
	.uleb128 0xf
	.long	0x32de
	.uleb128 0x29
	.byte	0x1
	.long	.LASF643
	.byte	0x4
	.byte	0x9
	.byte	0x1
	.long	0x33e
	.quad	.LFB62
	.quad	.LFE62
	.byte	0x1
	.byte	0x9c
	.long	0x3910
	.uleb128 0x24
	.string	"A"
	.byte	0x4
	.byte	0x9
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2a
	.long	.LASF637
	.byte	0x4
	.byte	0x9
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -172
	.uleb128 0x24
	.string	"is"
	.byte	0x4
	.byte	0x9
	.long	0xeb0
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x24
	.string	"ov"
	.byte	0x4
	.byte	0x9
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -188
	.uleb128 0x2b
	.string	"a"
	.byte	0x4
	.byte	0xb
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x2c
	.long	.LASF634
	.byte	0x4
	.byte	0xc
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2b
	.string	"row"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2b
	.string	"i"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2b
	.string	"j"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2b
	.string	"k"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2b
	.string	"l"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x2b
	.string	"m"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2b
	.string	"n"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -116
	.uleb128 0x2c
	.long	.LASF638
	.byte	0x4
	.byte	0xd
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2b
	.string	"isz"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -76
	.uleb128 0x2b
	.string	"val"
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2c
	.long	.LASF639
	.byte	0x4
	.byte	0xd
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2b
	.string	"idx"
	.byte	0x4
	.byte	0xe
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2c
	.long	.LASF640
	.byte	0x4
	.byte	0xf
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2b
	.string	"end"
	.byte	0x4
	.byte	0xf
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2b
	.string	"ai"
	.byte	0x4
	.byte	0xf
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2b
	.string	"aj"
	.byte	0x4
	.byte	0xf
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x2b
	.string	"bs"
	.byte	0x4
	.byte	0xf
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x2c
	.long	.LASF641
	.byte	0x4
	.byte	0xf
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2c
	.long	.LASF642
	.byte	0x4
	.byte	0x10
	.long	0x34f1
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x27
	.long	.LASF632
	.long	0x3926
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.24258
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x34e6
	.uleb128 0xd
	.long	0x126
	.long	0x3926
	.uleb128 0xe
	.long	0xe0
	.byte	0x1a
	.byte	0x0
	.uleb128 0xf
	.long	0x3916
	.uleb128 0x29
	.byte	0x1
	.long	.LASF644
	.byte	0x4
	.byte	0x4c
	.byte	0x1
	.long	0x33e
	.quad	.LFB63
	.quad	.LFE63
	.byte	0x1
	.byte	0x9c
	.long	0x3b3c
	.uleb128 0x24
	.string	"A"
	.byte	0x4
	.byte	0x4c
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2a
	.long	.LASF645
	.byte	0x4
	.byte	0x4c
	.long	0xd94
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x2a
	.long	.LASF646
	.byte	0x4
	.byte	0x4c
	.long	0xd94
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2a
	.long	.LASF647
	.byte	0x4
	.byte	0x4c
	.long	0x1b8c
	.byte	0x3
	.byte	0x91
	.sleb128 -236
	.uleb128 0x24
	.string	"B"
	.byte	0x4
	.byte	0x4c
	.long	0x2899
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x2b
	.string	"a"
	.byte	0x4
	.byte	0x4e
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2b
	.string	"c"
	.byte	0x4
	.byte	0x4e
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2c
	.long	.LASF634
	.byte	0x4
	.byte	0x4f
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2c
	.long	.LASF648
	.byte	0x4
	.byte	0x50
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2b
	.string	"i"
	.byte	0x4
	.byte	0x50
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -124
	.uleb128 0x2b
	.string	"k"
	.byte	0x4
	.byte	0x50
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2c
	.long	.LASF649
	.byte	0x4
	.byte	0x50
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -116
	.uleb128 0x2c
	.long	.LASF650
	.byte	0x4
	.byte	0x50
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x2c
	.long	.LASF651
	.byte	0x4
	.byte	0x50
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -108
	.uleb128 0x2c
	.long	.LASF652
	.byte	0x4
	.byte	0x50
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2b
	.string	"row"
	.byte	0x4
	.byte	0x51
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2c
	.long	.LASF653
	.byte	0x4
	.byte	0x51
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2c
	.long	.LASF654
	.byte	0x4
	.byte	0x51
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2c
	.long	.LASF655
	.byte	0x4
	.byte	0x51
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x2c
	.long	.LASF656
	.byte	0x4
	.byte	0x51
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2c
	.long	.LASF657
	.byte	0x4
	.byte	0x52
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2c
	.long	.LASF619
	.byte	0x4
	.byte	0x52
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2c
	.long	.LASF453
	.byte	0x4
	.byte	0x53
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -180
	.uleb128 0x2c
	.long	.LASF658
	.byte	0x4
	.byte	0x53
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2c
	.long	.LASF659
	.byte	0x4
	.byte	0x53
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2b
	.string	"bs"
	.byte	0x4
	.byte	0x53
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2b
	.string	"bs2"
	.byte	0x4
	.byte	0x53
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2b
	.string	"aj"
	.byte	0x4
	.byte	0x54
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2b
	.string	"ai"
	.byte	0x4
	.byte	0x54
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x2c
	.long	.LASF660
	.byte	0x4
	.byte	0x55
	.long	0x32ee
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x2b
	.string	"C"
	.byte	0x4
	.byte	0x56
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2c
	.long	.LASF661
	.byte	0x4
	.byte	0x57
	.long	0x3e0
	.byte	0x3
	.byte	0x91
	.sleb128 -196
	.uleb128 0x2c
	.long	.LASF662
	.byte	0x4
	.byte	0x57
	.long	0x3e0
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x27
	.long	.LASF632
	.long	0x3b3c
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.24627
	.byte	0x0
	.uleb128 0xf
	.long	0x6f5
	.uleb128 0x29
	.byte	0x1
	.long	.LASF663
	.byte	0x4
	.byte	0xa2
	.byte	0x1
	.long	0x33e
	.quad	.LFB64
	.quad	.LFE64
	.byte	0x1
	.byte	0x9c
	.long	0x3c78
	.uleb128 0x24
	.string	"A"
	.byte	0x4
	.byte	0xa2
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2a
	.long	.LASF645
	.byte	0x4
	.byte	0xa2
	.long	0xd94
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2a
	.long	.LASF646
	.byte	0x4
	.byte	0xa2
	.long	0xd94
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2a
	.long	.LASF647
	.byte	0x4
	.byte	0xa2
	.long	0x1b8c
	.byte	0x3
	.byte	0x91
	.sleb128 -140
	.uleb128 0x24
	.string	"B"
	.byte	0x4
	.byte	0xa2
	.long	0x2899
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2b
	.string	"a"
	.byte	0x4
	.byte	0xa4
	.long	0x3910
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2b
	.string	"is1"
	.byte	0x4
	.byte	0xa5
	.long	0xd94
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2b
	.string	"is2"
	.byte	0x4
	.byte	0xa5
	.long	0xd94
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2c
	.long	.LASF634
	.byte	0x4
	.byte	0xa6
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x2c
	.long	.LASF664
	.byte	0x4
	.byte	0xa7
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2c
	.long	.LASF665
	.byte	0x4
	.byte	0xa7
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2c
	.long	.LASF453
	.byte	0x4
	.byte	0xa7
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2c
	.long	.LASF658
	.byte	0x4
	.byte	0xa7
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2b
	.string	"i"
	.byte	0x4
	.byte	0xa7
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -44
	.uleb128 0x2b
	.string	"bs"
	.byte	0x4
	.byte	0xa7
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x2c
	.long	.LASF4
	.byte	0x4
	.byte	0xa7
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x2c
	.long	.LASF657
	.byte	0x4
	.byte	0xa8
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2c
	.long	.LASF619
	.byte	0x4
	.byte	0xa8
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x27
	.long	.LASF632
	.long	0x3c88
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.25037
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x3c88
	.uleb128 0xe
	.long	0xe0
	.byte	0x17
	.byte	0x0
	.uleb128 0xf
	.long	0x3c78
	.uleb128 0x29
	.byte	0x1
	.long	.LASF666
	.byte	0x4
	.byte	0xd0
	.byte	0x1
	.long	0x33e
	.quad	.LFB65
	.quad	.LFE65
	.byte	0x1
	.byte	0x9c
	.long	0x3d31
	.uleb128 0x24
	.string	"A"
	.byte	0x4
	.byte	0xd0
	.long	0x195a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x24
	.string	"n"
	.byte	0x4
	.byte	0xd0
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2a
	.long	.LASF657
	.byte	0x4
	.byte	0xd0
	.long	0x2aa4
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2a
	.long	.LASF619
	.byte	0x4
	.byte	0xd0
	.long	0x2aa4
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2a
	.long	.LASF647
	.byte	0x4
	.byte	0xd0
	.long	0x1b8c
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x24
	.string	"B"
	.byte	0x4
	.byte	0xd0
	.long	0x2aaf
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2c
	.long	.LASF634
	.byte	0x4
	.byte	0xd2
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x2b
	.string	"i"
	.byte	0x4
	.byte	0xd3
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x3d41
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.25328
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x3d41
	.uleb128 0xe
	.long	0xe0
	.byte	0x19
	.byte	0x0
	.uleb128 0xf
	.long	0x3d31
	.uleb128 0x29
	.byte	0x1
	.long	.LASF667
	.byte	0x4
	.byte	0xe7
	.byte	0x1
	.long	0x33e
	.quad	.LFB66
	.quad	.LFE66
	.byte	0x1
	.byte	0x9c
	.long	0x3ef2
	.uleb128 0x24
	.string	"A"
	.byte	0x4
	.byte	0xe7
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x24
	.string	"xx"
	.byte	0x4
	.byte	0xe7
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x24
	.string	"zz"
	.byte	0x4
	.byte	0xe7
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2b
	.string	"a"
	.byte	0x4
	.byte	0xe9
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2b
	.string	"z"
	.byte	0x4
	.byte	0xea
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2b
	.string	"sum"
	.byte	0x4
	.byte	0xea
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2b
	.string	"x"
	.byte	0x4
	.byte	0xeb
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2b
	.string	"v"
	.byte	0x4
	.byte	0xec
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2c
	.long	.LASF634
	.byte	0x4
	.byte	0xed
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -108
	.uleb128 0x2b
	.string	"mbs"
	.byte	0x4
	.byte	0xee
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2b
	.string	"i"
	.byte	0x4
	.byte	0xee
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2b
	.string	"n"
	.byte	0x4
	.byte	0xee
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2c
	.long	.LASF668
	.byte	0x4
	.byte	0xee
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2b
	.string	"idx"
	.byte	0x4
	.byte	0xef
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2b
	.string	"ii"
	.byte	0x4
	.byte	0xef
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2c
	.long	.LASF669
	.byte	0x4
	.byte	0xef
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2c
	.long	.LASF670
	.byte	0x4
	.byte	0xf0
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x27
	.long	.LASF632
	.long	0x3f0d
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.25437
	.uleb128 0x2d
	.quad	.LBB2
	.quad	.LBE2
	.long	0x3e9d
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x105
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x105
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.byte	0x0
	.uleb128 0x2d
	.quad	.LBB3
	.quad	.LBE3
	.long	0x3ed0
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x106
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x106
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB4
	.quad	.LBE4
	.uleb128 0x2e
	.string	"__i"
	.byte	0x4
	.value	0x108
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.byte	0x0
	.byte	0x0
	.uleb128 0x4
	.byte	0x8
	.long	0x3ef8
	.uleb128 0xf
	.long	0x3be
	.uleb128 0xd
	.long	0x126
	.long	0x3f0d
	.uleb128 0xe
	.long	0xe0
	.byte	0x11
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF672
	.byte	0x4
	.value	0x118
	.byte	0x1
	.long	0x33e
	.quad	.LFB67
	.quad	.LFE67
	.byte	0x1
	.byte	0x9c
	.long	0x4104
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x118
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x118
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x118
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x11a
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x11b
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x11b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x11b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x11b
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x11c
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x11c
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x11d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x11d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x11e
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x11f
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x120
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x120
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x120
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x120
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x120
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x120
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x120
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x120
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x121
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x4104
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.25630
	.uleb128 0x2d
	.quad	.LBB5
	.quad	.LBE5
	.long	0x40d4
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x137
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x137
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB6
	.quad	.LBE6
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x138
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x138
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF676
	.byte	0x4
	.value	0x14b
	.byte	0x1
	.long	0x33e
	.quad	.LFB68
	.quad	.LFE68
	.byte	0x1
	.byte	0x9c
	.long	0x431a
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x14b
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x14b
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x14b
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x14d
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x14e
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x14e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x14e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x14e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x14e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x14e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x14e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x14e
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x14f
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x14f
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x150
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x151
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x152
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x152
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x152
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x152
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x152
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x152
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x152
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x152
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x153
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x431a
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.25809
	.uleb128 0x2d
	.quad	.LBB7
	.quad	.LBE7
	.long	0x42ea
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x16e
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x16e
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB8
	.quad	.LBE8
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x16f
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x16f
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF678
	.byte	0x4
	.value	0x183
	.byte	0x1
	.long	0x33e
	.quad	.LFB69
	.quad	.LFE69
	.byte	0x1
	.byte	0x9c
	.long	0x454f
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x183
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x183
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x183
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x185
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x186
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x186
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x186
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x187
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x187
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x188
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x189
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x18a
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x18a
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x18a
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x18a
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x18a
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x18a
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x18a
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x18a
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x18b
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x454f
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.26011
	.uleb128 0x2d
	.quad	.LBB9
	.quad	.LBE9
	.long	0x451f
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x1a1
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x1a1
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB10
	.quad	.LBE10
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x1a2
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x1a2
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF680
	.byte	0x4
	.value	0x1b8
	.byte	0x1
	.long	0x33e
	.quad	.LFB70
	.quad	.LFE70
	.byte	0x1
	.byte	0x9c
	.long	0x47a2
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x1b8
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x1b8
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x1b8
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x1ba
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x1bb
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x1bb
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x1bb
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x1bc
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x1bc
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x1bd
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x1be
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x1bf
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x1bf
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x1bf
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x1c0
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x1c0
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x1c0
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x1c0
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x1c0
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x1c1
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x47a2
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.26244
	.uleb128 0x2d
	.quad	.LBB11
	.quad	.LBE11
	.long	0x4772
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x1d7
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x1d7
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB12
	.quad	.LBE12
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x1d8
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x1d8
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF682
	.byte	0x4
	.value	0x1f0
	.byte	0x1
	.long	0x33e
	.quad	.LFB71
	.quad	.LFE71
	.byte	0x1
	.byte	0x9c
	.long	0x4a15
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x1f0
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x1f0
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x1f0
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x1f2
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x1f3
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x1f3
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x1f3
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x1f3
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x1f3
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x1f3
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x1f3
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x1f4
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x1f4
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x1f5
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x1f5
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x1f5
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x1f5
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x1f5
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x6"
	.byte	0x4
	.value	0x1f5
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x1f5
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x1f6
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x1f7
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x1f8
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x1f8
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x1f8
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x1f8
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x1f8
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x1f8
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x1f8
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x1f8
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x1f9
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x4a15
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.26516
	.uleb128 0x2d
	.quad	.LBB13
	.quad	.LBE13
	.long	0x49e5
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x20f
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x20f
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB14
	.quad	.LBE14
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x210
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x210
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF684
	.byte	0x4
	.value	0x229
	.byte	0x1
	.long	0x33e
	.quad	.LFB72
	.quad	.LFE72
	.byte	0x1
	.byte	0x9c
	.long	0x4caa
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x229
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x229
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -304
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x229
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x22b
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x22c
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x22c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x22c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x22c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x22c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x22c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x22c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF685
	.byte	0x4
	.value	0x22c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x22d
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x22d
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x22e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x22e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x22e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x22e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x22e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x6"
	.byte	0x4
	.value	0x22e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x7"
	.byte	0x4
	.value	0x22e
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x22e
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x22f
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x230
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -116
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x231
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x231
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -108
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x231
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x231
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x231
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x231
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x231
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x231
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x232
	.long	0x3e0
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x27
	.long	.LASF632
	.long	0x4caa
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.26835
	.uleb128 0x2d
	.quad	.LBB15
	.quad	.LBE15
	.long	0x4c7a
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x248
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x248
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB16
	.quad	.LBE16
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x249
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x249
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF686
	.byte	0x4
	.value	0x266
	.byte	0x1
	.long	0x33e
	.quad	.LFB73
	.quad	.LFE73
	.byte	0x1
	.byte	0x9c
	.long	0x4f19
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x266
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x266
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -304
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x266
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x268
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x269
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF685
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF687
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF688
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF689
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF690
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x28
	.long	.LASF691
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x28
	.long	.LASF692
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x28
	.long	.LASF693
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x28
	.long	.LASF694
	.byte	0x4
	.value	0x269
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x26a
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x26a
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x26b
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x2e
	.string	"xv"
	.byte	0x4
	.value	0x26b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x26c
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x26d
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x26e
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"ij"
	.byte	0x4
	.value	0x26e
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x26e
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x26f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x26f
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x26f
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x26f
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x26f
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x26f
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x26f
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x270
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x4f19
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.27211
	.byte	0x0
	.uleb128 0xf
	.long	0x3c78
	.uleb128 0x30
	.byte	0x1
	.long	.LASF695
	.byte	0x4
	.value	0x2af
	.byte	0x1
	.long	0x33e
	.quad	.LFB74
	.quad	.LFE74
	.byte	0x1
	.byte	0x9c
	.long	0x51a7
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x2af
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x2af
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -320
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x2af
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -328
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x2b1
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x2b2
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x28
	.long	.LASF685
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF687
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF688
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF689
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF690
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF691
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF692
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x28
	.long	.LASF693
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x28
	.long	.LASF694
	.byte	0x4
	.value	0x2b2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x2b3
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -288
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x2b3
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x2b4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x2b4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x2b4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x2b4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x2b4
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x2b5
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x2b6
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x2b7
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ij"
	.byte	0x4
	.value	0x2b7
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x2b7
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x2b8
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x2b8
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x2b8
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x2b8
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x2b8
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x2b8
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x2b9
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x51a7
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.27446
	.byte	0x0
	.uleb128 0xf
	.long	0x3c78
	.uleb128 0x30
	.byte	0x1
	.long	.LASF696
	.byte	0x4
	.value	0x32e
	.byte	0x1
	.long	0x33e
	.quad	.LFB75
	.quad	.LFE75
	.byte	0x1
	.byte	0x9c
	.long	0x5471
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x32e
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -344
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x32e
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -352
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x32e
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -360
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x330
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x331
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -304
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -288
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x28
	.long	.LASF685
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x28
	.long	.LASF687
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x28
	.long	.LASF688
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x28
	.long	.LASF689
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x28
	.long	.LASF690
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF691
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF692
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF693
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF694
	.byte	0x4
	.value	0x331
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x332
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -320
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x332
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x6"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x7"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"x8"
	.byte	0x4
	.value	0x333
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x333
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -328
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x334
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x335
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x336
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ij"
	.byte	0x4
	.value	0x336
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x336
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x337
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x337
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x337
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x337
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x337
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x337
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x338
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x5471
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.28485
	.byte	0x0
	.uleb128 0xf
	.long	0x3c78
	.uleb128 0x30
	.byte	0x1
	.long	.LASF697
	.byte	0x4
	.value	0x38a
	.byte	0x1
	.long	0x33e
	.quad	.LFB76
	.quad	.LFE76
	.byte	0x1
	.byte	0x9c
	.long	0x57aa
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x38a
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -392
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x38a
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -400
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x38a
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -408
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x38c
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -368
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x38d
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -360
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -352
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -344
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -336
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -328
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -320
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x28
	.long	.LASF685
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -304
	.uleb128 0x28
	.long	.LASF687
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x28
	.long	.LASF688
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -288
	.uleb128 0x28
	.long	.LASF689
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x28
	.long	.LASF690
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x28
	.long	.LASF691
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x28
	.long	.LASF692
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x28
	.long	.LASF693
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x28
	.long	.LASF694
	.byte	0x4
	.value	0x38d
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x38e
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -376
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x38e
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"x6"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"x7"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x8"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x9"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x10"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x11"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x12"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x13"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x14"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"x15"
	.byte	0x4
	.value	0x38f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x38f
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -384
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x390
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x391
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x392
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ij"
	.byte	0x4
	.value	0x392
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x392
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x393
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x393
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x393
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x393
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x393
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x393
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x394
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x57aa
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.29563
	.byte	0x0
	.uleb128 0xf
	.long	0x3c78
	.uleb128 0x30
	.byte	0x1
	.long	.LASF698
	.byte	0x4
	.value	0x3d5
	.byte	0x1
	.long	0x33e
	.quad	.LFB77
	.quad	.LFE77
	.byte	0x1
	.byte	0x9c
	.long	0x59bc
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x3d5
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x3d5
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x3d5
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x3d7
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x3d8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x3d8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x3d8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x28
	.long	.LASF699
	.byte	0x4
	.value	0x3d8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF700
	.byte	0x4
	.value	0x3d8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x3d8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x3d9
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x3da
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x3db
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x3db
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x3db
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x3db
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"bs"
	.byte	0x4
	.value	0x3db
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x3db
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x3db
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"bs2"
	.byte	0x4
	.value	0x3db
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x28
	.long	.LASF658
	.byte	0x4
	.value	0x3dc
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x3dc
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x3dc
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x3dc
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x3dd
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x59bc
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.30632
	.uleb128 0x2f
	.quad	.LBB17
	.quad	.LBE17
	.uleb128 0x28
	.long	.LASF701
	.byte	0x4
	.value	0x3ff
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF702
	.byte	0x4
	.value	0x3ff
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF703
	.byte	0x4
	.value	0x3ff
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -188
	.uleb128 0x28
	.long	.LASF704
	.byte	0x4
	.value	0x3ff
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF705
	.byte	0x4
	.value	0x3ff
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -196
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x3efd
	.uleb128 0x30
	.byte	0x1
	.long	.LASF706
	.byte	0x4
	.value	0x40c
	.byte	0x1
	.long	0x33e
	.quad	.LFB78
	.quad	.LFE78
	.byte	0x1
	.byte	0x9c
	.long	0x5b9c
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x40c
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x40c
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x40c
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x40c
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x40e
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x40f
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"y"
	.byte	0x4
	.value	0x410
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x410
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"sum"
	.byte	0x4
	.value	0x410
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x411
	.long	0x3ef2
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x412
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x413
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -108
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x413
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x413
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x413
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x28
	.long	.LASF668
	.byte	0x4
	.value	0x413
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x414
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x414
	.long	0x15b7
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x415
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x27
	.long	.LASF632
	.long	0x5bac
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.30821
	.uleb128 0x2d
	.quad	.LBB18
	.quad	.LBE18
	.long	0x5b47
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x436
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x436
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.byte	0x0
	.uleb128 0x2d
	.quad	.LBB19
	.quad	.LBE19
	.long	0x5b7a
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x437
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x437
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB20
	.quad	.LBE20
	.uleb128 0x2e
	.string	"__i"
	.byte	0x4
	.value	0x438
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.byte	0x0
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x5bac
	.uleb128 0xe
	.long	0xe0
	.byte	0x14
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF707
	.byte	0x4
	.value	0x44c
	.byte	0x1
	.long	0x33e
	.quad	.LFB79
	.quad	.LFE79
	.byte	0x1
	.byte	0x9c
	.long	0x5dc1
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x44c
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x44c
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x44c
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x44c
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x44e
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x44f
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"y"
	.byte	0x4
	.value	0x44f
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x44f
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x44f
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x44f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x44f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x450
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x450
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF708
	.byte	0x4
	.value	0x450
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x450
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x451
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x452
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x453
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x453
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x453
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x453
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x453
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x453
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x453
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x454
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x5dc1
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.31048
	.uleb128 0x2d
	.quad	.LBB21
	.quad	.LBE21
	.long	0x5d91
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x478
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x478
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB22
	.quad	.LBE22
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x479
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x479
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF709
	.byte	0x4
	.value	0x490
	.byte	0x1
	.long	0x33e
	.quad	.LFB80
	.quad	.LFE80
	.byte	0x1
	.byte	0x9c
	.long	0x5ff5
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x490
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x490
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x490
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x490
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x492
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x493
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"y"
	.byte	0x4
	.value	0x493
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x493
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x493
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x493
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x493
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x493
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x493
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x493
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x493
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF708
	.byte	0x4
	.value	0x493
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x493
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x494
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x495
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x496
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x496
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x496
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x496
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x496
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x496
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x496
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x497
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x5ff5
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.31276
	.uleb128 0x2d
	.quad	.LBB23
	.quad	.LBE23
	.long	0x5fc5
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x4b8
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x4b8
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB24
	.quad	.LBE24
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x4b9
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x4b9
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF710
	.byte	0x4
	.value	0x4d1
	.byte	0x1
	.long	0x33e
	.quad	.LFB81
	.quad	.LFE81
	.byte	0x1
	.byte	0x9c
	.long	0x6248
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x4d1
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x4d1
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x4d1
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x4d1
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x4d3
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x4d4
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"y"
	.byte	0x4
	.value	0x4d4
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x4d4
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x4d4
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x4d4
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF708
	.byte	0x4
	.value	0x4d4
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x4d4
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x4d5
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x4d6
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x4d7
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x4d7
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x4d7
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x4d7
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x4d7
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x4d7
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x4d7
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x4d8
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x6248
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.31515
	.uleb128 0x2d
	.quad	.LBB25
	.quad	.LBE25
	.long	0x6218
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x4f9
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x4f9
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB26
	.quad	.LBE26
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x4fa
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x4fa
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF711
	.byte	0x4
	.value	0x514
	.byte	0x1
	.long	0x33e
	.quad	.LFB82
	.quad	.LFE82
	.byte	0x1
	.byte	0x9c
	.long	0x64ba
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x514
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x514
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x514
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x514
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -288
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x516
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x517
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"y"
	.byte	0x4
	.value	0x517
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x517
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x517
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x517
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF708
	.byte	0x4
	.value	0x518
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x518
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x519
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x51a
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x51b
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x51b
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x51b
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x51b
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x51b
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x51b
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x51b
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x51c
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x64ba
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.31786
	.uleb128 0x2d
	.quad	.LBB27
	.quad	.LBE27
	.long	0x648a
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x53d
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x53d
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB28
	.quad	.LBE28
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x53e
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x53e
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF712
	.byte	0x4
	.value	0x558
	.byte	0x1
	.long	0x33e
	.quad	.LFB83
	.quad	.LFE83
	.byte	0x1
	.byte	0x9c
	.long	0x674b
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x558
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x558
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -288
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x558
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x558
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -304
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x55a
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x55b
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x2e
	.string	"y"
	.byte	0x4
	.value	0x55b
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x55b
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x55b
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x55b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x55b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x55b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x55b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x55b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x55b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x55c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x55c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x55c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x55c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x55c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"x6"
	.byte	0x4
	.value	0x55c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF708
	.byte	0x4
	.value	0x55c
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x55c
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x55d
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x55e
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x55f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x55f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x55f
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x55f
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x55f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x55f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x55f
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x560
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x674b
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.32097
	.uleb128 0x2d
	.quad	.LBB29
	.quad	.LBE29
	.long	0x671b
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x581
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x581
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB30
	.quad	.LBE30
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x582
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x582
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF713
	.byte	0x4
	.value	0x59e
	.byte	0x1
	.long	0x33e
	.quad	.LFB84
	.quad	.LFE84
	.byte	0x1
	.byte	0x9c
	.long	0x69fd
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x59e
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x59e
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -320
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x59e
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -328
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x59e
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -336
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x5a0
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x5a1
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x2e
	.string	"y"
	.byte	0x4
	.value	0x5a1
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -264
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x5a1
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x5a1
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -248
	.uleb128 0x28
	.long	.LASF673
	.byte	0x4
	.value	0x5a1
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -240
	.uleb128 0x28
	.long	.LASF674
	.byte	0x4
	.value	0x5a1
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x28
	.long	.LASF677
	.byte	0x4
	.value	0x5a1
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x28
	.long	.LASF679
	.byte	0x4
	.value	0x5a1
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x28
	.long	.LASF681
	.byte	0x4
	.value	0x5a1
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x28
	.long	.LASF683
	.byte	0x4
	.value	0x5a1
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x28
	.long	.LASF685
	.byte	0x4
	.value	0x5a1
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x5a2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x5a2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x5a2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x5a2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x5a2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x6"
	.byte	0x4
	.value	0x5a2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"x7"
	.byte	0x4
	.value	0x5a2
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x28
	.long	.LASF708
	.byte	0x4
	.value	0x5a2
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -288
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x5a2
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x5a3
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x5a4
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -116
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x5a5
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x5a5
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -108
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x5a5
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x5a5
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x5a5
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x5a5
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x5a5
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x5a6
	.long	0x3e0
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x27
	.long	.LASF632
	.long	0x69fd
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.32456
	.uleb128 0x2d
	.quad	.LBB31
	.quad	.LBE31
	.long	0x69cd
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x5c7
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x5c7
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.byte	0x0
	.uleb128 0x2f
	.quad	.LBB32
	.quad	.LBE32
	.uleb128 0x2e
	.string	"_p"
	.byte	0x4
	.value	0x5c8
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF671
	.byte	0x4
	.value	0x5c8
	.long	0x333
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF714
	.byte	0x4
	.value	0x5e5
	.byte	0x1
	.long	0x33e
	.quad	.LFB85
	.quad	.LFE85
	.byte	0x1
	.byte	0x9c
	.long	0x6bff
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x5e5
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x5e5
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x5e5
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x5e5
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x5e7
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x5e8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x5e8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x5e8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x28
	.long	.LASF699
	.byte	0x4
	.value	0x5e8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x28
	.long	.LASF700
	.byte	0x4
	.value	0x5e8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF675
	.byte	0x4
	.value	0x5e8
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x5e9
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x5ea
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x5eb
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x5eb
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x5eb
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x5eb
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"bs"
	.byte	0x4
	.value	0x5eb
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x5eb
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x5eb
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"bs2"
	.byte	0x4
	.value	0x5eb
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x28
	.long	.LASF658
	.byte	0x4
	.value	0x5ec
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x5ec
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x5ec
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x5ed
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x6bff
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.32859
	.uleb128 0x2f
	.quad	.LBB33
	.quad	.LBE33
	.uleb128 0x28
	.long	.LASF701
	.byte	0x4
	.value	0x60f
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x28
	.long	.LASF703
	.byte	0x4
	.value	0x60f
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -180
	.uleb128 0x28
	.long	.LASF704
	.byte	0x4
	.value	0x60f
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x28
	.long	.LASF705
	.byte	0x4
	.value	0x60f
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -188
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x5b9c
	.uleb128 0x30
	.byte	0x1
	.long	.LASF715
	.byte	0x4
	.value	0x61e
	.byte	0x1
	.long	0x33e
	.quad	.LFB86
	.quad	.LFE86
	.byte	0x1
	.byte	0x9c
	.long	0x6c85
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x61e
	.long	0x195a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x61e
	.long	0xece
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x61e
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x28
	.long	.LASF716
	.byte	0x4
	.value	0x620
	.long	0x3a7
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x621
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x6c95
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.33035
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x6c95
	.uleb128 0xe
	.long	0xe0
	.byte	0x21
	.byte	0x0
	.uleb128 0xf
	.long	0x6c85
	.uleb128 0x30
	.byte	0x1
	.long	.LASF717
	.byte	0x4
	.value	0x62b
	.byte	0x1
	.long	0x33e
	.quad	.LFB87
	.quad	.LFE87
	.byte	0x1
	.byte	0x9c
	.long	0x6d1b
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x62b
	.long	0x195a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x62b
	.long	0xece
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x62b
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x28
	.long	.LASF716
	.byte	0x4
	.value	0x62d
	.long	0x3a7
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x62e
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x6d2b
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.33103
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x6d2b
	.uleb128 0xe
	.long	0xe0
	.byte	0x18
	.byte	0x0
	.uleb128 0xf
	.long	0x6d1b
	.uleb128 0x30
	.byte	0x1
	.long	.LASF718
	.byte	0x4
	.value	0x638
	.byte	0x1
	.long	0x33e
	.quad	.LFB88
	.quad	.LFE88
	.byte	0x1
	.byte	0x9c
	.long	0x6f75
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x638
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -296
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x638
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -304
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x638
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x638
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -320
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x63b
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"zb"
	.byte	0x4
	.value	0x63c
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x63c
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x63c
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x63c
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x63c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x63c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x63c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x63c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x63c
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x63d
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x63e
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -124
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x63f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x63f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -116
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x63f
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x63f
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x28
	.long	.LASF719
	.byte	0x4
	.value	0x63f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"bs"
	.byte	0x4
	.value	0x63f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x63f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x63f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"bs2"
	.byte	0x4
	.value	0x63f
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -76
	.uleb128 0x2e
	.string	"ib"
	.byte	0x4
	.value	0x63f
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x63f
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x28
	.long	.LASF720
	.byte	0x4
	.value	0x640
	.long	0x32b2
	.byte	0x3
	.byte	0x91
	.sleb128 -256
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x641
	.long	0x3e0
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x27
	.long	.LASF632
	.long	0x6f85
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.33195
	.uleb128 0x2f
	.quad	.LBB34
	.quad	.LBE34
	.uleb128 0x28
	.long	.LASF658
	.byte	0x4
	.value	0x6a7
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x6a7
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -44
	.uleb128 0x28
	.long	.LASF699
	.byte	0x4
	.value	0x6a8
	.long	0x9a3
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF700
	.byte	0x4
	.value	0x6a8
	.long	0x9a3
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF721
	.byte	0x4
	.value	0x6a8
	.long	0x9a3
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.byte	0x0
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x6f85
	.uleb128 0xe
	.long	0xe0
	.byte	0x24
	.byte	0x0
	.uleb128 0xf
	.long	0x6f75
	.uleb128 0x30
	.byte	0x1
	.long	.LASF722
	.byte	0x4
	.value	0x6cd
	.byte	0x1
	.long	0x33e
	.quad	.LFB89
	.quad	.LFE89
	.byte	0x1
	.byte	0x9c
	.long	0x7223
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x6cd
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -312
	.uleb128 0x26
	.string	"xx"
	.byte	0x4
	.value	0x6cd
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -320
	.uleb128 0x26
	.string	"yy"
	.byte	0x4
	.value	0x6cd
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -328
	.uleb128 0x26
	.string	"zz"
	.byte	0x4
	.value	0x6cd
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -336
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x6cf
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -216
	.uleb128 0x2e
	.string	"zb"
	.byte	0x4
	.value	0x6d0
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -208
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x6d0
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -224
	.uleb128 0x2e
	.string	"z"
	.byte	0x4
	.value	0x6d0
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -232
	.uleb128 0x2e
	.string	"xb"
	.byte	0x4
	.value	0x6d0
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -200
	.uleb128 0x2e
	.string	"x1"
	.byte	0x4
	.value	0x6d0
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -192
	.uleb128 0x2e
	.string	"x2"
	.byte	0x4
	.value	0x6d0
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"x3"
	.byte	0x4
	.value	0x6d0
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x2e
	.string	"x4"
	.byte	0x4
	.value	0x6d0
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x2e
	.string	"x5"
	.byte	0x4
	.value	0x6d0
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -160
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x6d1
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x6d2
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -140
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x6d3
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x6d3
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -132
	.uleb128 0x2e
	.string	"idx"
	.byte	0x4
	.value	0x6d3
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"ii"
	.byte	0x4
	.value	0x6d3
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x28
	.long	.LASF719
	.byte	0x4
	.value	0x6d3
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -108
	.uleb128 0x2e
	.string	"bs"
	.byte	0x4
	.value	0x6d3
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x6d3
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x6d3
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"bs2"
	.byte	0x4
	.value	0x6d3
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -92
	.uleb128 0x2e
	.string	"ib"
	.byte	0x4
	.value	0x6d3
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x28
	.long	.LASF669
	.byte	0x4
	.value	0x6d3
	.long	0x991
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x28
	.long	.LASF720
	.byte	0x4
	.value	0x6d4
	.long	0x32b2
	.byte	0x3
	.byte	0x91
	.sleb128 -272
	.uleb128 0x28
	.long	.LASF670
	.byte	0x4
	.value	0x6d5
	.long	0x3e0
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x27
	.long	.LASF632
	.long	0x7233
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.33932
	.uleb128 0x2f
	.quad	.LBB35
	.quad	.LBE35
	.uleb128 0x28
	.long	.LASF658
	.byte	0x4
	.value	0x73b
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x73b
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x28
	.long	.LASF699
	.byte	0x4
	.value	0x73c
	.long	0x9a3
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF700
	.byte	0x4
	.value	0x73c
	.long	0x9a3
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF721
	.byte	0x4
	.value	0x73c
	.long	0x9a3
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x2f
	.quad	.LBB36
	.quad	.LBE36
	.uleb128 0x28
	.long	.LASF701
	.byte	0x4
	.value	0x74b
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -280
	.uleb128 0x28
	.long	.LASF703
	.byte	0x4
	.value	0x74b
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -284
	.uleb128 0x28
	.long	.LASF704
	.byte	0x4
	.value	0x74b
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -288
	.uleb128 0x28
	.long	.LASF705
	.byte	0x4
	.value	0x74b
	.long	0x354
	.byte	0x3
	.byte	0x91
	.sleb128 -292
	.byte	0x0
	.byte	0x0
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x7233
	.uleb128 0xe
	.long	0xe0
	.byte	0x1b
	.byte	0x0
	.uleb128 0xf
	.long	0x7223
	.uleb128 0x30
	.byte	0x1
	.long	.LASF723
	.byte	0x4
	.value	0x760
	.byte	0x1
	.long	0x33e
	.quad	.LFB90
	.quad	.LFE90
	.byte	0x1
	.byte	0x9c
	.long	0x72e9
	.uleb128 0x26
	.string	"inA"
	.byte	0x4
	.value	0x760
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x31
	.long	.LASF469
	.byte	0x4
	.value	0x760
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x762
	.long	0x3910
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF724
	.byte	0x4
	.value	0x763
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x28
	.long	.LASF725
	.byte	0x4
	.value	0x764
	.long	0x3a7
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x765
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x2e
	.string	"one"
	.byte	0x4
	.value	0x766
	.long	0x354
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2e
	.string	"tnz"
	.byte	0x4
	.value	0x766
	.long	0x354
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x27
	.long	.LASF632
	.long	0x72f9
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.34647
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x72f9
	.uleb128 0xe
	.long	0xe0
	.byte	0x10
	.byte	0x0
	.uleb128 0xf
	.long	0x72e9
	.uleb128 0x30
	.byte	0x1
	.long	.LASF726
	.byte	0x4
	.value	0x770
	.byte	0x1
	.long	0x33e
	.quad	.LFB91
	.quad	.LFE91
	.byte	0x1
	.byte	0x9c
	.long	0x7433
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x770
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x31
	.long	.LASF79
	.byte	0x4
	.value	0x770
	.long	0xf98
	.byte	0x3
	.byte	0x91
	.sleb128 -124
	.uleb128 0x31
	.long	.LASF255
	.byte	0x4
	.value	0x770
	.long	0x824
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x772
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x773
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x774
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"sum"
	.byte	0x4
	.value	0x775
	.long	0x39c
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x776
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x776
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x776
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2e
	.string	"bs"
	.byte	0x4
	.value	0x776
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"nz"
	.byte	0x4
	.value	0x776
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x2e
	.string	"bs2"
	.byte	0x4
	.value	0x776
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x2e
	.string	"k1"
	.byte	0x4
	.value	0x776
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -44
	.uleb128 0x27
	.long	.LASF632
	.long	0x7433
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.34724
	.uleb128 0x2f
	.quad	.LBB37
	.quad	.LBE37
	.uleb128 0x2e
	.string	"tmp"
	.byte	0x4
	.value	0x783
	.long	0x824
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF727
	.byte	0x4
	.value	0x784
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.long	0x32de
	.uleb128 0x30
	.byte	0x1
	.long	.LASF728
	.byte	0x4
	.value	0x7ab
	.byte	0x1
	.long	0x33e
	.quad	.LFB92
	.quad	.LFE92
	.byte	0x1
	.byte	0x9c
	.long	0x74c6
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x7ab
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x26
	.string	"B"
	.byte	0x4
	.value	0x7ab
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x26
	.string	"flg"
	.byte	0x4
	.value	0x7ab
	.long	0x28e5
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x7ad
	.long	0x3910
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"b"
	.byte	0x4
	.value	0x7ad
	.long	0x3910
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x7ae
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x27
	.long	.LASF632
	.long	0x74c6
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.34937
	.byte	0x0
	.uleb128 0xf
	.long	0x72e9
	.uleb128 0x30
	.byte	0x1
	.long	.LASF729
	.byte	0x4
	.value	0x7ca
	.byte	0x1
	.long	0x33e
	.quad	.LFB93
	.quad	.LFE93
	.byte	0x1
	.byte	0x9c
	.long	0x7609
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x7ca
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x26
	.string	"v"
	.byte	0x4
	.value	0x7ca
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x7cc
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x7cd
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -84
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -76
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -100
	.uleb128 0x2e
	.string	"row"
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"bs"
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2e
	.string	"ai"
	.byte	0x4
	.value	0x7ce
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"aj"
	.byte	0x4
	.value	0x7ce
	.long	0x991
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x28
	.long	.LASF730
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x2e
	.string	"bs2"
	.byte	0x4
	.value	0x7ce
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -44
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x7cf
	.long	0x9a3
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x28
	.long	.LASF716
	.byte	0x4
	.value	0x7cf
	.long	0x3a7
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x2e
	.string	"aa"
	.byte	0x4
	.value	0x7d0
	.long	0x32ee
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF731
	.byte	0x4
	.value	0x7d0
	.long	0x32ee
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.uleb128 0x27
	.long	.LASF632
	.long	0x7619
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.35139
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x7619
	.uleb128 0xe
	.long	0xe0
	.byte	0x16
	.byte	0x0
	.uleb128 0xf
	.long	0x7609
	.uleb128 0x30
	.byte	0x1
	.long	.LASF732
	.byte	0x4
	.value	0x7ef
	.byte	0x1
	.long	0x33e
	.quad	.LFB94
	.quad	.LFE94
	.byte	0x1
	.byte	0x9c
	.long	0x77db
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x7ef
	.long	0x195a
	.byte	0x3
	.byte	0x91
	.sleb128 -168
	.uleb128 0x26
	.string	"ll"
	.byte	0x4
	.value	0x7ef
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -176
	.uleb128 0x26
	.string	"rr"
	.byte	0x4
	.value	0x7ef
	.long	0xece
	.byte	0x3
	.byte	0x91
	.sleb128 -184
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x7f1
	.long	0x3910
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x2e
	.string	"l"
	.byte	0x4
	.value	0x7f2
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x2e
	.string	"r"
	.byte	0x4
	.value	0x7f2
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x2e
	.string	"li"
	.byte	0x4
	.value	0x7f2
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x2e
	.string	"ri"
	.byte	0x4
	.value	0x7f2
	.long	0x150e
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x2e
	.string	"x"
	.byte	0x4
	.value	0x7f3
	.long	0x3a7
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x2e
	.string	"aa"
	.byte	0x4
	.value	0x7f4
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x2e
	.string	"v"
	.byte	0x4
	.value	0x7f4
	.long	0x32ee
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x7f5
	.long	0x33e
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x2e
	.string	"i"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -76
	.uleb128 0x2e
	.string	"j"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x2e
	.string	"k"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -68
	.uleb128 0x2e
	.string	"lm"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -148
	.uleb128 0x2e
	.string	"rn"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x2e
	.string	"M"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x2e
	.string	"m"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -60
	.uleb128 0x2e
	.string	"n"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"mbs"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -52
	.uleb128 0x2e
	.string	"tmp"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x2e
	.string	"bs"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -44
	.uleb128 0x2e
	.string	"bs2"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x2e
	.string	"iai"
	.byte	0x4
	.value	0x7f6
	.long	0x36a
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x2e
	.string	"ai"
	.byte	0x4
	.value	0x7f7
	.long	0x15b7
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x2e
	.string	"aj"
	.byte	0x4
	.value	0x7f7
	.long	0x15b7
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.uleb128 0x27
	.long	.LASF632
	.long	0x77db
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.35288
	.byte	0x0
	.uleb128 0xf
	.long	0x6d1b
	.uleb128 0x30
	.byte	0x1
	.long	.LASF733
	.byte	0x4
	.value	0x82e
	.byte	0x1
	.long	0x33e
	.quad	.LFB95
	.quad	.LFE95
	.byte	0x1
	.byte	0x9c
	.long	0x7851
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x82e
	.long	0x195a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x31
	.long	.LASF661
	.byte	0x4
	.value	0x82e
	.long	0x1d63
	.byte	0x2
	.byte	0x91
	.sleb128 -44
	.uleb128 0x31
	.long	.LASF347
	.byte	0x4
	.value	0x82e
	.long	0x28bf
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x830
	.long	0x3910
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.uleb128 0x27
	.long	.LASF632
	.long	0x7861
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.35495
	.byte	0x0
	.uleb128 0xd
	.long	0x126
	.long	0x7861
	.uleb128 0xe
	.long	0xe0
	.byte	0x12
	.byte	0x0
	.uleb128 0xf
	.long	0x7851
	.uleb128 0x30
	.byte	0x1
	.long	.LASF734
	.byte	0x4
	.value	0x849
	.byte	0x1
	.long	0x33e
	.quad	.LFB96
	.quad	.LFE96
	.byte	0x1
	.byte	0x9c
	.long	0x78c8
	.uleb128 0x26
	.string	"A"
	.byte	0x4
	.value	0x849
	.long	0x195a
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x2e
	.string	"a"
	.byte	0x4
	.value	0x84b
	.long	0x3910
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x28
	.long	.LASF634
	.byte	0x4
	.value	0x84c
	.long	0x33e
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x27
	.long	.LASF632
	.long	0x78c8
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	__func__.35575
	.byte	0x0
	.uleb128 0xf
	.long	0x7609
	.uleb128 0x32
	.long	.LASF735
	.byte	0x15
	.byte	0xa7
	.long	0x307
	.byte	0x1
	.byte	0x1
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x7903
	.uleb128 0x1d
	.long	0xd5
	.uleb128 0x1d
	.long	0x2d
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0xc8
	.byte	0x0
	.uleb128 0x33
	.long	.LASF736
	.byte	0x2
	.value	0x45b
	.long	0x7911
	.byte	0x1
	.byte	0x1
	.uleb128 0x4
	.byte	0x8
	.long	0x78da
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x793b
	.uleb128 0x1d
	.long	0x4b
	.uleb128 0x1d
	.long	0x2d
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x1d
	.long	0x333
	.byte	0x0
	.uleb128 0x33
	.long	.LASF737
	.byte	0x2
	.value	0x45c
	.long	0x7949
	.byte	0x1
	.byte	0x1
	.uleb128 0x4
	.byte	0x8
	.long	0x7917
	.uleb128 0x33
	.long	.LASF738
	.byte	0xd
	.value	0x18d
	.long	0x795d
	.byte	0x1
	.byte	0x1
	.uleb128 0x4
	.byte	0x8
	.long	0x800
	.uleb128 0x1c
	.byte	0x1
	.long	0x33e
	.long	0x7974
	.uleb128 0x1d
	.long	0x333
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x33
	.long	.LASF739
	.byte	0x2
	.value	0x638
	.long	0x7982
	.byte	0x1
	.byte	0x1
	.uleb128 0x4
	.byte	0x8
	.long	0x7963
	.uleb128 0x32
	.long	.LASF740
	.byte	0xe
	.byte	0x21
	.long	0x3b2
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF741
	.byte	0xe
	.byte	0xa2
	.long	0xcca
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF742
	.byte	0x14
	.byte	0x20
	.long	0x126
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF743
	.byte	0x14
	.byte	0x21
	.long	0x126
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF744
	.byte	0x14
	.byte	0x22
	.long	0x36a
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF735
	.byte	0x15
	.byte	0xa7
	.long	0x307
	.byte	0x1
	.byte	0x1
	.uleb128 0x33
	.long	.LASF736
	.byte	0x2
	.value	0x45b
	.long	0x7911
	.byte	0x1
	.byte	0x1
	.uleb128 0x33
	.long	.LASF737
	.byte	0x2
	.value	0x45c
	.long	0x7949
	.byte	0x1
	.byte	0x1
	.uleb128 0x33
	.long	.LASF738
	.byte	0xd
	.value	0x18d
	.long	0x795d
	.byte	0x1
	.byte	0x1
	.uleb128 0x33
	.long	.LASF739
	.byte	0x2
	.value	0x638
	.long	0x7982
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF740
	.byte	0xe
	.byte	0x21
	.long	0x3b2
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF741
	.byte	0xe
	.byte	0xa2
	.long	0xcca
	.byte	0x1
	.byte	0x1
	.uleb128 0x35
	.long	.LASF745
	.byte	0x16
	.byte	0x2e
	.long	0x2d
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	mtNexts
	.uleb128 0x35
	.long	.LASF746
	.byte	0x16
	.byte	0x2f
	.long	0x6e3
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	mtNexti
	.uleb128 0xd
	.long	0x6e3
	.long	0x7a65
	.uleb128 0x36
	.long	0xe0
	.value	0x26f
	.byte	0x0
	.uleb128 0x35
	.long	.LASF747
	.byte	0x16
	.byte	0x30
	.long	0x7a54
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	s_seeds
	.uleb128 0x32
	.long	.LASF742
	.byte	0x14
	.byte	0x20
	.long	0x126
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF743
	.byte	0x14
	.byte	0x21
	.long	0x126
	.byte	0x1
	.byte	0x1
	.uleb128 0x32
	.long	.LASF744
	.byte	0x14
	.byte	0x22
	.long	0x36a
	.byte	0x1
	.byte	0x1
	.byte	0x0
	.section	.debug_abbrev
	.uleb128 0x1
	.uleb128 0x11
	.byte	0x1
	.uleb128 0x25
	.uleb128 0xe
	.uleb128 0x13
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x1b
	.uleb128 0xe
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x10
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x2
	.uleb128 0x24
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0x8
	.byte	0x0
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x16
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x4
	.uleb128 0xf
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x5
	.uleb128 0xf
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.byte	0x0
	.byte	0x0
	.uleb128 0x6
	.uleb128 0x16
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x7
	.uleb128 0x24
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.byte	0x0
	.byte	0x0
	.uleb128 0x8
	.uleb128 0x13
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x9
	.uleb128 0xd
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xd
	.byte	0x0
	.byte	0x0
	.uleb128 0xa
	.uleb128 0x16
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.byte	0x0
	.byte	0x0
	.uleb128 0xb
	.uleb128 0x13
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0xc
	.uleb128 0xd
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xd
	.byte	0x0
	.byte	0x0
	.uleb128 0xd
	.uleb128 0x1
	.byte	0x1
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0xe
	.uleb128 0x21
	.byte	0x0
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2f
	.uleb128 0xb
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.uleb128 0x26
	.byte	0x0
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x10
	.uleb128 0x4
	.byte	0x1
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x11
	.uleb128 0x28
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x1c
	.uleb128 0xd
	.byte	0x0
	.byte	0x0
	.uleb128 0x12
	.uleb128 0x4
	.byte	0x1
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x13
	.uleb128 0x13
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0x5
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x14
	.uleb128 0xd
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xd
	.byte	0x0
	.byte	0x0
	.uleb128 0x15
	.uleb128 0x13
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3c
	.uleb128 0xc
	.byte	0x0
	.byte	0x0
	.uleb128 0x16
	.uleb128 0x26
	.byte	0x0
	.byte	0x0
	.byte	0x0
	.uleb128 0x17
	.uleb128 0x4
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x18
	.uleb128 0x15
	.byte	0x0
	.uleb128 0x27
	.uleb128 0xc
	.byte	0x0
	.byte	0x0
	.uleb128 0x19
	.uleb128 0x15
	.byte	0x0
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1a
	.uleb128 0x13
	.byte	0x1
	.uleb128 0xb
	.uleb128 0x5
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1b
	.uleb128 0x13
	.byte	0x1
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1c
	.uleb128 0x15
	.byte	0x1
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1d
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1e
	.uleb128 0x16
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1f
	.uleb128 0xd
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xd
	.byte	0x0
	.byte	0x0
	.uleb128 0x20
	.uleb128 0x13
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0x5
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x21
	.uleb128 0x13
	.byte	0x1
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x22
	.uleb128 0x13
	.byte	0x1
	.uleb128 0xb
	.uleb128 0x5
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x23
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0xa
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x24
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x25
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0xa
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x26
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x27
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x34
	.uleb128 0xc
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x28
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x29
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0xa
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x2a
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x2b
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x2c
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x2d
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x2e
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x2f
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.byte	0x0
	.byte	0x0
	.uleb128 0x30
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0xa
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x31
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x32
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3c
	.uleb128 0xc
	.byte	0x0
	.byte	0x0
	.uleb128 0x33
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3c
	.uleb128 0xc
	.byte	0x0
	.byte	0x0
	.uleb128 0x34
	.uleb128 0x18
	.byte	0x0
	.byte	0x0
	.byte	0x0
	.uleb128 0x35
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x36
	.uleb128 0x21
	.byte	0x0
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2f
	.uleb128 0x5
	.byte	0x0
	.byte	0x0
	.byte	0x0
	.section	.debug_pubnames,"",@progbits
	.long	0x3cb
	.value	0x2
	.long	.Ldebug_info0
	.long	0x7aa3
	.long	0x377b
	.string	"MatIncreaseOverlap_SeqBAIJ"
	.long	0x392b
	.string	"MatGetSubMatrix_SeqBAIJ_Private"
	.long	0x3b41
	.string	"MatGetSubMatrix_SeqBAIJ"
	.long	0x3c8d
	.string	"MatGetSubMatrices_SeqBAIJ"
	.long	0x3d46
	.string	"MatMult_SeqBAIJ_1"
	.long	0x3f12
	.string	"MatMult_SeqBAIJ_2"
	.long	0x4109
	.string	"MatMult_SeqBAIJ_3"
	.long	0x431f
	.string	"MatMult_SeqBAIJ_4"
	.long	0x4554
	.string	"MatMult_SeqBAIJ_5"
	.long	0x47a7
	.string	"MatMult_SeqBAIJ_6"
	.long	0x4a1a
	.string	"MatMult_SeqBAIJ_7"
	.long	0x4caf
	.string	"MatMult_SeqBAIJ_15_ver1"
	.long	0x4f1e
	.string	"MatMult_SeqBAIJ_15_ver2"
	.long	0x51ac
	.string	"MatMult_SeqBAIJ_15_ver3"
	.long	0x5476
	.string	"MatMult_SeqBAIJ_15_ver4"
	.long	0x57af
	.string	"MatMult_SeqBAIJ_N"
	.long	0x59c1
	.string	"MatMultAdd_SeqBAIJ_1"
	.long	0x5bb1
	.string	"MatMultAdd_SeqBAIJ_2"
	.long	0x5dc6
	.string	"MatMultAdd_SeqBAIJ_3"
	.long	0x5ffa
	.string	"MatMultAdd_SeqBAIJ_4"
	.long	0x624d
	.string	"MatMultAdd_SeqBAIJ_5"
	.long	0x64bf
	.string	"MatMultAdd_SeqBAIJ_6"
	.long	0x6750
	.string	"MatMultAdd_SeqBAIJ_7"
	.long	0x6a02
	.string	"MatMultAdd_SeqBAIJ_N"
	.long	0x6c04
	.string	"MatMultHermitianTranspose_SeqBAIJ"
	.long	0x6c9a
	.string	"MatMultTranspose_SeqBAIJ"
	.long	0x6d30
	.string	"MatMultHermitianTransposeAdd_SeqBAIJ"
	.long	0x6f8a
	.string	"MatMultTransposeAdd_SeqBAIJ"
	.long	0x7238
	.string	"MatScale_SeqBAIJ"
	.long	0x72fe
	.string	"MatNorm_SeqBAIJ"
	.long	0x7438
	.string	"MatEqual_SeqBAIJ"
	.long	0x74cb
	.string	"MatGetDiagonal_SeqBAIJ"
	.long	0x761e
	.string	"MatDiagonalScale_SeqBAIJ"
	.long	0x77e0
	.string	"MatGetInfo_SeqBAIJ"
	.long	0x7866
	.string	"MatZeroEntries_SeqBAIJ"
	.long	0x7a28
	.string	"mtNexts"
	.long	0x7a3e
	.string	"mtNexti"
	.long	0x7a65
	.string	"s_seeds"
	.long	0x0
	.section	.debug_pubtypes,"",@progbits
	.long	0x627
	.value	0x2
	.long	.Ldebug_info0
	.long	0x7aa3
	.long	0x34
	.string	"MPI_Comm"
	.long	0x4d
	.string	"MPI_Request"
	.long	0x67
	.string	"MPI_Status"
	.long	0xb6
	.string	"MPI_Status"
	.long	0xd5
	.string	"size_t"
	.long	0x10a
	.string	"__off_t"
	.long	0x115
	.string	"__off64_t"
	.long	0x2c9
	.string	"_IO_lock_t"
	.long	0x2d0
	.string	"_IO_marker"
	.long	0x12d
	.string	"_IO_FILE"
	.long	0x33e
	.string	"PetscErrorCode"
	.long	0x349
	.string	"PetscClassId"
	.long	0x354
	.string	"PetscBLASInt"
	.long	0x35f
	.string	"PetscMPIInt"
	.long	0x36a
	.string	"PetscInt"
	.long	0x38a
	.string	"PetscPrecision"
	.long	0x39c
	.string	"PetscReal"
	.long	0x3a7
	.string	"PetscScalar"
	.long	0x3b2
	.string	"PetscLogDouble"
	.long	0x3be
	.string	"MatScalar"
	.long	0x3e0
	.string	"PetscBool"
	.long	0x408
	.string	"PetscObject"
	.long	0x6b4
	.string	"PetscFList"
	.long	0x6cc
	.string	"PetscViewer"
	.long	0x6e3
	.string	"uint"
	.long	0x705
	.string	"uintptr_t"
	.long	0x717
	.string	"H5F_mem_t"
	.long	0x800
	.string	"PetscStack"
	.long	0x80c
	.string	"PetscOList"
	.long	0x980
	.string	"PetscOps"
	.long	0x41a
	.string	"_p_PetscObject"
	.long	0xa00
	.string	"_p_PetscObject"
	.long	0xa0b
	.string	"PetscIntStack"
	.long	0xa43
	.string	"PetscClassRegInfo"
	.long	0xa92
	.string	"PetscClassPerfInfo"
	.long	0xa9d
	.string	"PetscClassRegLog"
	.long	0xaae
	.string	"_n_PetscClassRegLog"
	.long	0xae5
	.string	"PetscClassPerfLog"
	.long	0xaf6
	.string	"_n_PetscClassPerfLog"
	.long	0xb4e
	.string	"PetscEventRegInfo"
	.long	0xbd9
	.string	"PetscEventPerfInfo"
	.long	0xbe4
	.string	"PetscEventRegLog"
	.long	0xbf5
	.string	"_n_PetscEventRegLog"
	.long	0xc2c
	.string	"PetscEventPerfLog"
	.long	0xc3d
	.string	"_n_PetscEventPerfLog"
	.long	0xc74
	.string	"_PetscStageInfo"
	.long	0xcbf
	.string	"PetscStageInfo"
	.long	0xcca
	.string	"PetscStageLog"
	.long	0xcdb
	.string	"_n_PetscStageLog"
	.long	0xd42
	.string	"PetscRandom"
	.long	0xd88
	.string	"InsertMode"
	.long	0xd94
	.string	"IS"
	.long	0xdaa
	.string	"_p_ISLocalToGlobalMapping"
	.long	0xe10
	.string	"ISLocalToGlobalMapping"
	.long	0xe36
	.string	"ISColoringType"
	.long	0xe41
	.string	"ISColoringValue"
	.long	0xe4c
	.string	"_n_ISColoring"
	.long	0xebc
	.string	"ISColoring"
	.long	0xece
	.string	"Vec"
	.long	0xf98
	.string	"NormType"
	.long	0xfb9
	.string	"VecOption"
	.long	0xfc5
	.string	"PetscLayout"
	.long	0xfd6
	.string	"_n_PetscLayout"
	.long	0x105c
	.string	"_VecOps"
	.long	0x18e7
	.string	"VecStash"
	.long	0x1913
	.string	"PetscCUSPFlag"
	.long	0xedf
	.string	"_p_Vec"
	.long	0x1924
	.string	"_n_Vecs"
	.long	0x1948
	.string	"Vecs"
	.long	0x195a
	.string	"Mat"
	.long	0x1b66
	.string	"MatFactorType"
	.long	0x1b8c
	.string	"MatReuse"
	.long	0x1bb8
	.string	"MatStructure"
	.long	0x1bd9
	.string	"MatAssemblyType"
	.long	0x1c79
	.string	"MatOption"
	.long	0x1ca1
	.string	"MatDuplicateOption"
	.long	0x1d3b
	.string	"MatInfo"
	.long	0x1d63
	.string	"MatInfoType"
	.long	0x1e0a
	.string	"MatFactorInfo"
	.long	0x1e5e
	.string	"MatSORType"
	.long	0x1e6a
	.string	"MatFDColoring"
	.long	0x1ff7
	.string	"MatNullSpace"
	.long	0x2099
	.string	"_MatOps"
	.long	0x303b
	.string	"PetscMatStashSpace"
	.long	0x304c
	.string	"_MatStashSpace"
	.long	0x3204
	.string	"MatStash"
	.long	0x325d
	.string	"MatStencilInfo"
	.long	0x32b2
	.string	"Mat_CompressedRow"
	.long	0x196b
	.string	"_p_Mat"
	.long	0x1e7c
	.string	"_p_MatFDColoring"
	.long	0x2009
	.string	"_p_MatNullSpace"
	.long	0x34e6
	.string	"Mat_SeqBAIJ"
	.long	0x34f1
	.string	"PetscBT"
	.long	0x0
	.section	.debug_aranges,"",@progbits
	.long	0x2c
	.value	0x2
	.long	.Ldebug_info0
	.byte	0x8
	.byte	0x0
	.value	0x0
	.value	0x0
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.quad	0x0
	.quad	0x0
	.section	.debug_str,"MS",@progbits,1
.LASF323:
	.string	"rindices"
.LASF50:
	.string	"_unused2"
.LASF36:
	.string	"_fileno"
.LASF453:
	.string	"nrows"
.LASF632:
	.string	"__func__"
.LASF474:
	.string	"restorerow"
.LASF165:
	.string	"PetscClassPerfInfo"
.LASF207:
	.string	"INSERT_ALL_VALUES"
.LASF392:
	.string	"MAT_KEEP_NONZERO_PATTERN"
.LASF133:
	.string	"H5FD_MEM_BTREE"
.LASF468:
	.string	"vecs"
.LASF606:
	.string	"ilen"
.LASF9:
	.string	"MPI_Status"
.LASF516:
	.string	"permute"
.LASF322:
	.string	"sindices"
.LASF635:
	.string	"VecGetArray"
.LASF739:
	.string	"PetscErrorPrintf"
.LASF558:
	.string	"getrowuppertriangular"
.LASF451:
	.string	"ncolumns"
.LASF378:
	.string	"SAME_PRECONDITIONER"
.LASF375:
	.string	"DIFFERENT_NONZERO_PATTERN"
.LASF60:
	.string	"PETSC_PRECISION_SINGLE"
.LASF69:
	.string	"PETSC_TRUE"
.LASF120:
	.string	"optionctx"
.LASF297:
	.string	"getvalues"
.LASF210:
	.string	"_p_IS"
.LASF511:
	.string	"getcolumnij"
.LASF300:
	.string	"create"
.LASF41:
	.string	"_shortbuf"
.LASF136:
	.string	"H5FD_MEM_LHEAP"
.LASF510:
	.string	"restorerowij"
.LASF281:
	.string	"dot_local"
.LASF213:
	.string	"globalstart"
.LASF675:
	.string	"zarray"
.LASF64:
	.string	"PetscReal"
.LASF605:
	.string	"imax"
.LASF530:
	.string	"multconstrained"
.LASF655:
	.string	"tcol"
.LASF423:
	.string	"MAT_GLOBAL_MAX"
.LASF182:
	.string	"_n_PetscEventRegLog"
.LASF586:
	.string	"total_space_size"
.LASF102:
	.string	"real_idmax"
.LASF503:
	.string	"getsubmatrices"
.LASF364:
	.string	"MAT_FACTOR_NONE"
.LASF116:
	.string	"python_destroy"
.LASF90:
	.string	"parentid"
.LASF22:
	.string	"_flags"
.LASF550:
	.string	"ptapsymbolic_seqaij"
.LASF367:
	.string	"MAT_FACTOR_ILU"
.LASF584:
	.string	"next"
.LASF269:
	.string	"setvalues"
.LASF541:
	.string	"matmult"
.LASF188:
	.string	"_PetscStageInfo"
.LASF301:
	.string	"stridegather"
.LASF18:
	.string	"__off_t"
.LASF561:
	.string	"getredundantmatrix"
.LASF112:
	.string	"scalarcomposeddata"
.LASF746:
	.string	"mtNexti"
.LASF710:
	.string	"MatMultAdd_SeqBAIJ_4"
.LASF471:
	.string	"rmctx"
.LASF745:
	.string	"mtNexts"
.LASF419:
	.string	"fill_ratio_needed"
.LASF537:
	.string	"ishermitian"
.LASF611:
	.string	"xtoyB"
.LASF593:
	.string	"MatStash"
.LASF645:
	.string	"isrow"
.LASF42:
	.string	"_lock"
.LASF98:
	.string	"intcomposedstate"
.LASF137:
	.string	"H5FD_MEM_OHDR"
.LASF106:
	.string	"realcomposeddata"
.LASF232:
	.string	"spptr"
.LASF316:
	.string	"recv_waits"
.LASF437:
	.string	"SOR_FORWARD_SWEEP"
.LASF482:
	.string	"solvetransposeadd"
.LASF231:
	.string	"valid_GPU_array"
.LASF175:
	.string	"visible"
.LASF637:
	.string	"is_max"
.LASF238:
	.string	"NormType"
.LASF517:
	.string	"getsubmatrix"
.LASF689:
	.string	"sum10"
.LASF690:
	.string	"sum11"
.LASF691:
	.string	"sum12"
.LASF692:
	.string	"sum13"
.LASF693:
	.string	"sum14"
.LASF694:
	.string	"sum15"
.LASF127:
	.string	"uint"
.LASF266:
	.string	"axpbypcz"
.LASF544:
	.string	"ptap"
.LASF327:
	.string	"donotstash"
.LASF283:
	.string	"norm_local"
.LASF177:
	.string	"numMessages"
.LASF240:
	.string	"VEC_IGNORE_NEGATIVE_INDICES"
.LASF147:
	.string	"PetscStack"
.LASF506:
	.string	"diagonalset"
.LASF77:
	.string	"bops"
.LASF377:
	.string	"SAME_NONZERO_PATTERN"
.LASF702:
	.string	"_zero"
.LASF446:
	.string	"SOR_APPLY_LOWER"
.LASF526:
	.string	"setcoloring"
.LASF533:
	.string	"mults"
.LASF625:
	.string	"idiagvalid"
.LASF351:
	.string	"stencil"
.LASF315:
	.string	"send_waits"
.LASF235:
	.string	"NORM_FROBENIUS"
.LASF638:
	.string	"nidx"
.LASF273:
	.string	"getsize"
.LASF472:
	.string	"_MatOps"
.LASF543:
	.string	"matmultnumeric"
.LASF248:
	.string	"mapping"
.LASF88:
	.string	"type_name"
.LASF190:
	.string	"perfInfo"
.LASF723:
	.string	"MatScale_SeqBAIJ"
.LASF670:
	.string	"usecprow"
.LASF246:
	.string	"range"
.LASF215:
	.string	"globals"
.LASF28:
	.string	"_IO_write_end"
.LASF410:
	.string	"MatDuplicateOption"
.LASF668:
	.string	"nonzerorow"
.LASF75:
	.string	"_p_PetscObject"
.LASF237:
	.string	"NORM_1_AND_2"
.LASF360:
	.string	"nooffprocentries"
.LASF464:
	.string	"ftn_func_cntx"
.LASF490:
	.string	"zeroentries"
.LASF473:
	.string	"getrow"
.LASF559:
	.string	"restorerowuppertriangular"
.LASF338:
	.string	"Vecs"
.LASF416:
	.string	"assemblies"
.LASF344:
	.string	"was_assembled"
.LASF104:
	.string	"realcomposedstate"
.LASF256:
	.string	"tdot"
.LASF110:
	.string	"scalarcomposedstate"
.LASF229:
	.string	"bstash"
.LASF341:
	.string	"cmap"
.LASF450:
	.string	"ncolors"
.LASF293:
	.string	"maxpointwisedivide"
.LASF539:
	.string	"setvaluesblockedlocal"
.LASF193:
	.string	"PetscStageInfo"
.LASF609:
	.string	"ignorezeroentries"
.LASF411:
	.string	"block_size"
.LASF393:
	.string	"MAT_IGNORE_ZERO_ENTRIES"
.LASF601:
	.string	"nonew"
.LASF72:
	.string	"PETSC_OWN_POINTER"
.LASF665:
	.string	"iary"
.LASF65:
	.string	"PetscScalar"
.LASF436:
	.string	"MatFactorInfo"
.LASF629:
	.string	"PetscMemcpy"
.LASF74:
	.string	"PetscObject"
.LASF365:
	.string	"MAT_FACTOR_LU"
.LASF217:
	.string	"IS_COLORING_GLOBAL"
.LASF470:
	.string	"remove"
.LASF93:
	.string	"tablevel"
.LASF582:
	.string	"PetscMatStashSpace"
.LASF607:
	.string	"free_imax_ilen"
.LASF737:
	.string	"PetscTrFree"
.LASF577:
	.string	"dummy4"
.LASF553:
	.string	"ptapnumeric_mpiaij"
.LASF253:
	.string	"destroyvecs"
.LASF716:
	.string	"zero"
.LASF219:
	.string	"ISColoringType"
.LASF79:
	.string	"type"
.LASF571:
	.string	"multhermitiantranspose"
.LASF68:
	.string	"PETSC_FALSE"
.LASF476:
	.string	"multadd"
.LASF180:
	.string	"PetscEventPerfInfo"
.LASF452:
	.string	"columns"
.LASF557:
	.string	"imaginarypart"
.LASF730:
	.string	"ambs"
.LASF205:
	.string	"ADD_VALUES"
.LASF534:
	.string	"solves"
.LASF295:
	.string	"pointwisemaxabs"
.LASF604:
	.string	"maxnz"
.LASF494:
	.string	"choleskyfactorsymbolic"
.LASF271:
	.string	"assemblyend"
.LASF718:
	.string	"MatMultHermitianTransposeAdd_SeqBAIJ"
.LASF35:
	.string	"_chain"
.LASF372:
	.string	"MAT_REUSE_MATRIX"
.LASF82:
	.string	"refct"
.LASF347:
	.string	"info"
.LASF528:
	.string	"setvaluesadifor"
.LASF13:
	.string	"unsigned char"
.LASF575:
	.string	"getcolumnnorms"
.LASF438:
	.string	"SOR_BACKWARD_SWEEP"
.LASF647:
	.string	"scall"
.LASF535:
	.string	"getinertia"
.LASF330:
	.string	"bowners"
.LASF358:
	.string	"spd_set"
.LASF751:
	.string	"_IO_lock_t"
.LASF704:
	.string	"_bbs"
.LASF63:
	.string	"float"
.LASF5:
	.string	"cancelled"
.LASF202:
	.string	"_p_PetscRandom"
.LASF396:
	.string	"MAT_SYMMETRY_ETERNAL"
.LASF616:
	.string	"diag"
.LASF119:
	.string	"optiondestroy"
.LASF263:
	.string	"maxpy"
.LASF479:
	.string	"solve"
.LASF546:
	.string	"ptapnumeric"
.LASF218:
	.string	"IS_COLORING_GHOSTED"
.LASF509:
	.string	"getrowij"
.LASF94:
	.string	"amem"
.LASF146:
	.string	"currentsize"
.LASF59:
	.string	"PetscInt"
.LASF735:
	.string	"stderr"
.LASF380:
	.string	"MAT_FLUSH_ASSEMBLY"
.LASF245:
	.string	"rend"
.LASF326:
	.string	"nprocessed"
.LASF673:
	.string	"sum1"
.LASF674:
	.string	"sum2"
.LASF677:
	.string	"sum3"
.LASF679:
	.string	"sum4"
.LASF681:
	.string	"sum5"
.LASF683:
	.string	"sum6"
.LASF603:
	.string	"singlemalloc"
.LASF687:
	.string	"sum8"
.LASF688:
	.string	"sum9"
.LASF408:
	.string	"MAT_COPY_VALUES"
.LASF153:
	.string	"compose"
.LASF502:
	.string	"iccfactor"
.LASF61:
	.string	"PETSC_PRECISION_DOUBLE"
.LASF560:
	.string	"matsolve"
.LASF131:
	.string	"H5FD_MEM_DEFAULT"
.LASF251:
	.string	"duplicate"
.LASF743:
	.string	"_BT_c"
.LASF307:
	.string	"umax"
.LASF183:
	.string	"numEvents"
.LASF260:
	.string	"swap"
.LASF439:
	.string	"SOR_SYMMETRIC_SWEEP"
.LASF523:
	.string	"getrowmaxabs"
.LASF676:
	.string	"MatMult_SeqBAIJ_3"
.LASF201:
	.string	"PetscRandom"
.LASF163:
	.string	"destructions"
.LASF105:
	.string	"realstarcomposedstate"
.LASF345:
	.string	"num_ass"
.LASF524:
	.string	"getrowminabs"
.LASF627:
	.string	"PetscBT"
.LASF748:
	.string	"GNU C 4.4.6 20110731 (Red Hat 4.4.6-3)"
.LASF198:
	.string	"stack"
.LASF66:
	.string	"PetscLogDouble"
.LASF184:
	.string	"maxEvents"
.LASF361:
	.string	"nooffproczerorows"
.LASF477:
	.string	"multtranspose"
.LASF596:
	.string	"MatStencilInfo"
.LASF27:
	.string	"_IO_write_ptr"
.LASF542:
	.string	"matmultsymbolic"
.LASF495:
	.string	"choleskyfactornumeric"
.LASF7:
	.string	"MPI_TAG"
.LASF157:
	.string	"publish"
.LASF648:
	.string	"smap"
.LASF619:
	.string	"icol"
.LASF406:
	.string	"MatOption"
.LASF335:
	.string	"PETSC_CUSP_BOTH"
.LASF415:
	.string	"memory"
.LASF99:
	.string	"intstarcomposedstate"
.LASF204:
	.string	"INSERT_VALUES"
.LASF623:
	.string	"sbaijMat"
.LASF580:
	.string	"setstencil"
.LASF114:
	.string	"fortran_func_pointers"
.LASF166:
	.string	"PetscClassRegLog"
.LASF160:
	.string	"_n_PetscIntStack"
.LASF521:
	.string	"unscalesystem"
.LASF145:
	.string	"line"
.LASF456:
	.string	"error_rel"
.LASF284:
	.string	"mdot_local"
.LASF447:
	.string	"MatSORType"
.LASF384:
	.string	"MAT_NEW_NONZERO_LOCATIONS"
.LASF299:
	.string	"shift"
.LASF311:
	.string	"size"
.LASF749:
	.string	"baij2.c"
.LASF429:
	.string	"dtcount"
.LASF595:
	.string	"starts"
.LASF149:
	.string	"_n_PetscOList"
.LASF115:
	.string	"python_context"
.LASF57:
	.string	"PetscBLASInt"
.LASF254:
	.string	"mdot"
.LASF321:
	.string	"rvalues"
.LASF0:
	.string	"MPI_Comm"
.LASF196:
	.string	"numStages"
.LASF11:
	.string	"size_t"
.LASF140:
	.string	"PETSC_ERROR_REPEAT"
.LASF548:
	.string	"matmulttransposesymbolic"
.LASF108:
	.string	"scalar_idmax"
.LASF417:
	.string	"mallocs"
.LASF152:
	.string	"destroy"
.LASF742:
	.string	"_BT_mask"
.LASF569:
	.string	"multdiagonalblock"
.LASF270:
	.string	"assemblybegin"
.LASF654:
	.string	"mat_j"
.LASF589:
	.string	"space"
.LASF85:
	.string	"class_name"
.LASF250:
	.string	"_VecOps"
.LASF407:
	.string	"MAT_DO_NOT_COPY_VALUES"
.LASF402:
	.string	"MAT_SPD"
.LASF113:
	.string	"scalarstarcomposeddata"
.LASF618:
	.string	"solve_work"
.LASF31:
	.string	"_IO_save_base"
.LASF56:
	.string	"PetscClassId"
.LASF223:
	.string	"ctype"
.LASF418:
	.string	"fill_ratio_given"
.LASF97:
	.string	"intstar_idmax"
.LASF728:
	.string	"MatEqual_SeqBAIJ"
.LASF404:
	.string	"MAT_NO_OFF_PROC_ZERO_ROWS"
.LASF454:
	.string	"rows"
.LASF551:
	.string	"ptapnumeric_seqaij"
.LASF413:
	.string	"nz_used"
.LASF504:
	.string	"increaseoverlap"
.LASF214:
	.string	"globalend"
.LASF491:
	.string	"zerorows"
.LASF568:
	.string	"restorelocalsubmatrix"
.LASF514:
	.string	"coloringpatch"
.LASF289:
	.string	"setlocaltoglobalmapping"
.LASF143:
	.string	"file"
.LASF741:
	.string	"_stageLog"
.LASF45:
	.string	"__pad2"
.LASF430:
	.string	"fill"
.LASF349:
	.string	"nearnullsp"
.LASF208:
	.string	"ADD_ALL_VALUES"
.LASF141:
	.string	"PETSC_ERROR_IN_CXX"
.LASF457:
	.string	"umin"
.LASF83:
	.string	"qlist"
.LASF164:
	.string	"descMem"
.LASF394:
	.string	"MAT_USE_INODES"
.LASF278:
	.string	"setvaluesblocked"
.LASF168:
	.string	"numClasses"
.LASF731:
	.string	"aa_j"
.LASF515:
	.string	"setunfactored"
.LASF722:
	.string	"MatMultTransposeAdd_SeqBAIJ"
.LASF414:
	.string	"nz_unneeded"
.LASF401:
	.string	"MAT_UNUSED_NONZERO_LOCATION_ERR"
.LASF602:
	.string	"nounused"
.LASF52:
	.string	"_next"
.LASF203:
	.string	"NOT_SET_VALUES"
.LASF513:
	.string	"fdcoloringcreate"
.LASF463:
	.string	"ftn_func_pointer"
.LASF151:
	.string	"view"
.LASF549:
	.string	"matmulttransposenumeric"
.LASF590:
	.string	"flg_v"
.LASF597:
	.string	"check"
.LASF628:
	.string	"PetscAbsScalar"
.LASF594:
	.string	"dims"
.LASF527:
	.string	"setvaluesadic"
.LASF650:
	.string	"kend"
.LASF257:
	.string	"mtdot"
.LASF486:
	.string	"getinfo"
.LASF306:
	.string	"nmax"
.LASF328:
	.string	"ignorenegidx"
.LASF130:
	.string	"H5FD_MEM_NOLIST"
.LASF374:
	.string	"MatReuse"
.LASF126:
	.string	"_p_PetscViewer"
.LASF614:
	.string	"free_a"
.LASF135:
	.string	"H5FD_MEM_GHEAP"
.LASF617:
	.string	"free_diag"
.LASF158:
	.string	"PetscOps"
.LASF682:
	.string	"MatMult_SeqBAIJ_6"
.LASF620:
	.string	"mult_work"
.LASF103:
	.string	"realstar_idmax"
.LASF738:
	.string	"petscstack"
.LASF531:
	.string	"multtransposeconstrained"
.LASF252:
	.string	"duplicatevecs"
.LASF159:
	.string	"PetscIntStack"
.LASF556:
	.string	"realpart"
.LASF334:
	.string	"PETSC_CUSP_CPU"
.LASF80:
	.string	"flops"
.LASF669:
	.string	"ridx"
.LASF379:
	.string	"MatStructure"
.LASF91:
	.string	"name"
.LASF525:
	.string	"convert"
.LASF53:
	.string	"_sbuf"
.LASF33:
	.string	"_IO_save_end"
.LASF443:
	.string	"SOR_ZERO_INITIAL_GUESS"
.LASF124:
	.string	"PetscViewer"
.LASF261:
	.string	"axpy"
.LASF661:
	.string	"flag"
.LASF505:
	.string	"getrowmax"
.LASF564:
	.string	"missingdiagonal"
.LASF715:
	.string	"MatMultHermitianTranspose_SeqBAIJ"
.LASF274:
	.string	"getlocalsize"
.LASF412:
	.string	"nz_allocated"
.LASF247:
	.string	"refcnt"
.LASF125:
	.string	"_n_PetscFList"
.LASF67:
	.string	"MatScalar"
.LASF636:
	.string	"VecRestoreArray"
.LASF244:
	.string	"rstart"
.LASF747:
	.string	"s_seeds"
.LASF6:
	.string	"MPI_SOURCE"
.LASF399:
	.string	"MAT_ERROR_LOWER_TRIANGULAR"
.LASF70:
	.string	"PetscBool"
.LASF729:
	.string	"MatGetDiagonal_SeqBAIJ"
.LASF441:
	.string	"SOR_LOCAL_BACKWARD_SWEEP"
.LASF14:
	.string	"short unsigned int"
.LASF16:
	.string	"signed char"
.LASF640:
	.string	"start"
.LASF195:
	.string	"_n_PetscStageLog"
.LASF455:
	.string	"columnsforrow"
.LASF325:
	.string	"nprocs"
.LASF428:
	.string	"dtcol"
.LASF78:
	.string	"comm"
.LASF4:
	.string	"count"
.LASF107:
	.string	"realstarcomposeddata"
.LASF422:
	.string	"MAT_LOCAL"
.LASF211:
	.string	"_p_ISLocalToGlobalMapping"
.LASF228:
	.string	"stash"
.LASF194:
	.string	"PetscStageLog"
.LASF275:
	.string	"restorearray"
.LASF613:
	.string	"free_ij"
.LASF320:
	.string	"svalues"
.LASF522:
	.string	"zerorowslocal"
.LASF305:
	.string	"restoresubvector"
.LASF19:
	.string	"__off64_t"
.LASF287:
	.string	"reciprocal"
.LASF25:
	.string	"_IO_read_base"
.LASF43:
	.string	"_offset"
.LASF425:
	.string	"MatInfoType"
.LASF302:
	.string	"stridescatter"
.LASF353:
	.string	"hermitian"
.LASF95:
	.string	"state"
.LASF30:
	.string	"_IO_buf_end"
.LASF55:
	.string	"PetscErrorCode"
.LASF615:
	.string	"compressedrow"
.LASF598:
	.string	"rindex"
.LASF1:
	.string	"MPI_Request"
.LASF660:
	.string	"mat_a"
.LASF227:
	.string	"array_gotten"
.LASF500:
	.string	"backwardsolve"
.LASF49:
	.string	"_mode"
.LASF262:
	.string	"axpby"
.LASF336:
	.string	"PetscCUSPFlag"
.LASF71:
	.string	"PETSC_COPY_VALUES"
.LASF282:
	.string	"tdot_local"
.LASF26:
	.string	"_IO_write_base"
.LASF405:
	.string	"NUM_MAT_OPTIONS"
.LASF142:
	.string	"function"
.LASF578:
	.string	"getsubmatricesparallel"
.LASF658:
	.string	"ncols"
.LASF461:
	.string	"currentcolor"
.LASF220:
	.string	"ISColoringValue"
.LASF386:
	.string	"MAT_STRUCTURALLY_SYMMETRIC"
.LASF721:
	.string	"xtmp"
.LASF359:
	.string	"symmetric_eternal"
.LASF234:
	.string	"NORM_2"
.LASF259:
	.string	"copy"
.LASF421:
	.string	"MatInfo"
.LASF352:
	.string	"symmetric"
.LASF81:
	.string	"time"
.LASF264:
	.string	"aypx"
.LASF308:
	.string	"oldnmax"
.LASF641:
	.string	"nidx2"
.LASF460:
	.string	"vscale"
.LASF431:
	.string	"levels"
.LASF272:
	.string	"getarray"
.LASF2:
	.string	"long int"
.LASF148:
	.string	"PetscOList"
.LASF579:
	.string	"setvaluesbatch"
.LASF685:
	.string	"sum7"
.LASF51:
	.string	"_IO_marker"
.LASF520:
	.string	"scalesystem"
.LASF501:
	.string	"ilufactor"
.LASF630:
	.string	"PetscMemzero"
.LASF701:
	.string	"_one"
.LASF329:
	.string	"insertmode"
.LASF381:
	.string	"MAT_FINAL_ASSEMBLY"
.LASF304:
	.string	"getsubvector"
.LASF545:
	.string	"ptapsymbolic"
.LASF646:
	.string	"iscol"
.LASF610:
	.string	"xtoy"
.LASF132:
	.string	"H5FD_MEM_SUPER"
.LASF331:
	.string	"VecStash"
.LASF562:
	.string	"getrowmin"
.LASF563:
	.string	"getcolumnvector"
.LASF488:
	.string	"getdiagonal"
.LASF483:
	.string	"lufactor"
.LASF591:
	.string	"reproduce"
.LASF371:
	.string	"MAT_INITIAL_MATRIX"
.LASF565:
	.string	"getseqnonzerostructure"
.LASF129:
	.string	"uintptr_t"
.LASF448:
	.string	"MatFDColoring"
.LASF176:
	.string	"depth"
.LASF622:
	.string	"saved_values"
.LASF356:
	.string	"hermitian_set"
.LASF12:
	.string	"long unsigned int"
.LASF744:
	.string	"_BT_idx"
.LASF376:
	.string	"SUBSET_NONZERO_PATTERN"
.LASF671:
	.string	"_end"
.LASF426:
	.string	"diagonal_fill"
.LASF538:
	.string	"isstructurallysymmetric"
.LASF294:
	.string	"pointwisemax"
.LASF599:
	.string	"Mat_CompressedRow"
.LASF651:
	.string	"oldcols"
.LASF686:
	.string	"MatMult_SeqBAIJ_15_ver1"
.LASF695:
	.string	"MatMult_SeqBAIJ_15_ver2"
.LASF696:
	.string	"MatMult_SeqBAIJ_15_ver3"
.LASF697:
	.string	"MatMult_SeqBAIJ_15_ver4"
.LASF309:
	.string	"reallocs"
.LASF409:
	.string	"MAT_SHARE_NONZERO_PATTERN"
.LASF191:
	.string	"eventLog"
.LASF171:
	.string	"PetscClassPerfLog"
.LASF216:
	.string	"ISLocalToGlobalMapping"
.LASF400:
	.string	"MAT_GETROW_UPPERTRIANGULAR"
.LASF20:
	.string	"char"
.LASF391:
	.string	"MAT_USE_HASH_TABLE"
.LASF362:
	.string	"valid_GPU_matrix"
.LASF621:
	.string	"sor_work"
.LASF343:
	.string	"assembled"
.LASF181:
	.string	"PetscEventRegLog"
.LASF122:
	.string	"optionsprinted"
.LASF186:
	.string	"PetscEventPerfLog"
.LASF111:
	.string	"scalarstarcomposedstate"
.LASF642:
	.string	"table"
.LASF280:
	.string	"replacearray"
.LASF58:
	.string	"PetscMPIInt"
.LASF29:
	.string	"_IO_buf_base"
.LASF187:
	.string	"_n_PetscEventPerfLog"
.LASF667:
	.string	"MatMult_SeqBAIJ_1"
.LASF672:
	.string	"MatMult_SeqBAIJ_2"
.LASF498:
	.string	"iccfactorsymbolic"
.LASF678:
	.string	"MatMult_SeqBAIJ_4"
.LASF680:
	.string	"MatMult_SeqBAIJ_5"
.LASF286:
	.string	"load"
.LASF684:
	.string	"MatMult_SeqBAIJ_7"
.LASF121:
	.string	"precision"
.LASF496:
	.string	"setuppreallocation"
.LASF206:
	.string	"MAX_VALUES"
.LASF178:
	.string	"messageLength"
.LASF390:
	.string	"MAT_NEW_NONZERO_ALLOCATION_ERR"
.LASF698:
	.string	"MatMult_SeqBAIJ_N"
.LASF118:
	.string	"optionhandler"
.LASF24:
	.string	"_IO_read_end"
.LASF317:
	.string	"send_status"
.LASF101:
	.string	"intstarcomposeddata"
.LASF167:
	.string	"_n_PetscClassRegLog"
.LASF719:
	.string	"rval"
.LASF382:
	.string	"MatAssemblyType"
.LASF21:
	.string	"_IO_FILE"
.LASF134:
	.string	"H5FD_MEM_DRAW"
.LASF489:
	.string	"diagonalscale"
.LASF109:
	.string	"scalarstar_idmax"
.LASF699:
	.string	"work"
.LASF663:
	.string	"MatGetSubMatrix_SeqBAIJ"
.LASF87:
	.string	"mansec"
.LASF139:
	.string	"PETSC_ERROR_INITIAL"
.LASF313:
	.string	"tag1"
.LASF314:
	.string	"tag2"
.LASF397:
	.string	"MAT_CHECK_COMPRESSED_ROW"
.LASF585:
	.string	"space_head"
.LASF332:
	.string	"PETSC_CUSP_UNALLOCATED"
.LASF466:
	.string	"_p_MatNullSpace"
.LASF459:
	.string	"vscaleforrow"
.LASF449:
	.string	"_p_MatFDColoring"
.LASF239:
	.string	"VEC_IGNORE_OFF_PROC_ENTRIES"
.LASF100:
	.string	"intcomposeddata"
.LASF276:
	.string	"setrandom"
.LASF497:
	.string	"ilufactorsymbolic"
.LASF403:
	.string	"MAT_NO_OFF_PROC_ENTRIES"
.LASF44:
	.string	"__pad1"
.LASF46:
	.string	"__pad3"
.LASF47:
	.string	"__pad4"
.LASF48:
	.string	"__pad5"
.LASF267:
	.string	"pointwisemult"
.LASF209:
	.string	"InsertMode"
.LASF292:
	.string	"setfromoptions"
.LASF445:
	.string	"SOR_APPLY_UPPER"
.LASF420:
	.string	"factor_mallocs"
.LASF653:
	.string	"mat_i"
.LASF34:
	.string	"_markers"
.LASF54:
	.string	"_pos"
.LASF708:
	.string	"yarray"
.LASF389:
	.string	"MAT_NEW_NONZERO_LOCATION_ERR"
.LASF138:
	.string	"H5FD_MEM_NTYPES"
.LASF62:
	.string	"PetscPrecision"
.LASF612:
	.string	"XtoY"
.LASF626:
	.string	"Mat_SeqBAIJ"
.LASF643:
	.string	"MatIncreaseOverlap_SeqBAIJ"
.LASF566:
	.string	"getghosts"
.LASF720:
	.string	"cprow"
.LASF155:
	.string	"composefunction"
.LASF96:
	.string	"int_idmax"
.LASF529:
	.string	"fdcoloringapply"
.LASF659:
	.string	"ssmap"
.LASF192:
	.string	"classLog"
.LASF10:
	.string	"double"
.LASF657:
	.string	"irow"
.LASF485:
	.string	"transpose"
.LASF624:
	.string	"idiag"
.LASF434:
	.string	"shifttype"
.LASF734:
	.string	"MatZeroEntries_SeqBAIJ"
.LASF512:
	.string	"restorecolumnij"
.LASF84:
	.string	"olist"
.LASF150:
	.string	"getcomm"
.LASF258:
	.string	"scale"
.LASF154:
	.string	"query"
.LASF633:
	.string	"VecRestoreArrayRead"
.LASF592:
	.string	"reproduce_count"
.LASF226:
	.string	"data"
.LASF76:
	.string	"classid"
.LASF705:
	.string	"_bncols"
.LASF123:
	.string	"PetscFList"
.LASF435:
	.string	"shiftamount"
.LASF296:
	.string	"pointwisemin"
.LASF222:
	.string	"colors"
.LASF197:
	.string	"maxStages"
.LASF117:
	.string	"noptionhandler"
.LASF644:
	.string	"MatGetSubMatrix_SeqBAIJ_Private"
.LASF241:
	.string	"VecOption"
.LASF363:
	.string	"solvertype"
.LASF732:
	.string	"MatDiagonalScale_SeqBAIJ"
.LASF161:
	.string	"PetscClassRegInfo"
.LASF540:
	.string	"getvecs"
.LASF574:
	.string	"findnonzerorows"
.LASF290:
	.string	"setvalueslocal"
.LASF703:
	.string	"_ione"
.LASF666:
	.string	"MatGetSubMatrices_SeqBAIJ"
.LASF173:
	.string	"PetscEventRegInfo"
.LASF717:
	.string	"MatMultTranspose_SeqBAIJ"
.LASF733:
	.string	"MatGetInfo_SeqBAIJ"
.LASF265:
	.string	"waxpy"
.LASF634:
	.string	"ierr"
.LASF225:
	.string	"_p_Vec"
.LASF480:
	.string	"solveadd"
.LASF340:
	.string	"rmap"
.LASF355:
	.string	"symmetric_set"
.LASF324:
	.string	"rmax"
.LASF576:
	.string	"invertblockdiagonal"
.LASF128:
	.string	"long long unsigned int"
.LASF89:
	.string	"parent"
.LASF481:
	.string	"solvetranspose"
.LASF467:
	.string	"has_cnst"
.LASF170:
	.string	"classInfo"
.LASF39:
	.string	"_cur_column"
.LASF740:
	.string	"_TotalFlops"
.LASF583:
	.string	"_MatStashSpace"
.LASF156:
	.string	"queryfunction"
.LASF656:
	.string	"mat_ilen"
.LASF427:
	.string	"usedt"
.LASF398:
	.string	"MAT_IGNORE_LOWER_TRIANGULAR"
.LASF337:
	.string	"_n_Vecs"
.LASF440:
	.string	"SOR_LOCAL_FORWARD_SWEEP"
.LASF444:
	.string	"SOR_EISENSTAT"
.LASF649:
	.string	"kstart"
.LASF277:
	.string	"setoption"
.LASF268:
	.string	"pointwisedivide"
.LASF354:
	.string	"structurally_symmetric"
.LASF319:
	.string	"nrecvs"
.LASF484:
	.string	"choleskyfactor"
.LASF726:
	.string	"MatNorm_SeqBAIJ"
.LASF32:
	.string	"_IO_backup_base"
.LASF387:
	.string	"MAT_NEW_DIAGONALS"
.LASF567:
	.string	"getlocalsubmatrix"
.LASF233:
	.string	"NORM_1"
.LASF23:
	.string	"_IO_read_ptr"
.LASF169:
	.string	"maxClasses"
.LASF552:
	.string	"ptapsymbolic_mpiaij"
.LASF185:
	.string	"eventInfo"
.LASF310:
	.string	"array"
.LASF706:
	.string	"MatMultAdd_SeqBAIJ_1"
.LASF707:
	.string	"MatMultAdd_SeqBAIJ_2"
.LASF709:
	.string	"MatMultAdd_SeqBAIJ_3"
.LASF442:
	.string	"SOR_LOCAL_SYMMETRIC_SWEEP"
.LASF711:
	.string	"MatMultAdd_SeqBAIJ_5"
.LASF712:
	.string	"MatMultAdd_SeqBAIJ_6"
.LASF713:
	.string	"MatMultAdd_SeqBAIJ_7"
.LASF144:
	.string	"directory"
.LASF588:
	.string	"local_remaining"
.LASF162:
	.string	"creations"
.LASF714:
	.string	"MatMultAdd_SeqBAIJ_N"
.LASF224:
	.string	"ISColoring"
.LASF493:
	.string	"lufactornumeric"
.LASF664:
	.string	"vary"
.LASF570:
	.string	"hermitiantranspose"
.LASF333:
	.string	"PETSC_CUSP_GPU"
.LASF312:
	.string	"rank"
.LASF279:
	.string	"placearray"
.LASF554:
	.string	"setsizes"
.LASF350:
	.string	"preallocated"
.LASF38:
	.string	"_old_offset"
.LASF547:
	.string	"matmulttranspose"
.LASF631:
	.string	"VecGetArrayRead"
.LASF469:
	.string	"alpha"
.LASF478:
	.string	"multtransposeadd"
.LASF179:
	.string	"numReductions"
.LASF532:
	.string	"findzerodiagonals"
.LASF339:
	.string	"_p_Mat"
.LASF291:
	.string	"resetarray"
.LASF388:
	.string	"MAT_IGNORE_OFF_PROC_ENTRIES"
.LASF3:
	.string	"long long int"
.LASF499:
	.string	"forwardsolve"
.LASF221:
	.string	"_n_ISColoring"
.LASF492:
	.string	"lufactorsymbolic"
.LASF346:
	.string	"same_nonzero"
.LASF366:
	.string	"MAT_FACTOR_CHOLESKY"
.LASF37:
	.string	"_flags2"
.LASF370:
	.string	"MatFactorType"
.LASF92:
	.string	"prefix"
.LASF395:
	.string	"MAT_HERMITIAN"
.LASF383:
	.string	"MAT_ROW_ORIENTED"
.LASF303:
	.string	"dotnorm2"
.LASF700:
	.string	"workt"
.LASF369:
	.string	"MAT_FACTOR_ILUDT"
.LASF200:
	.string	"stageInfo"
.LASF662:
	.string	"sorted"
.LASF368:
	.string	"MAT_FACTOR_ICC"
.LASF507:
	.string	"zerorowscolumns"
.LASF725:
	.string	"oalpha"
.LASF8:
	.string	"MPI_ERROR"
.LASF243:
	.string	"_n_PetscLayout"
.LASF189:
	.string	"used"
.LASF724:
	.string	"totalnz"
.LASF249:
	.string	"bmapping"
.LASF536:
	.string	"issymmetric"
.LASF752:
	.string	"H5F_mem_t"
.LASF639:
	.string	"ival"
.LASF465:
	.string	"MatNullSpace"
.LASF487:
	.string	"equal"
.LASF236:
	.string	"NORM_INFINITY"
.LASF581:
	.string	"setgrid"
.LASF573:
	.string	"getmultiprocblock"
.LASF230:
	.string	"petscnative"
.LASF750:
	.string	"/home/dpnkarthik/petsc-rnet/src/mat/impls/baij/seq"
.LASF357:
	.string	"structurally_symmetric_set"
.LASF255:
	.string	"norm"
.LASF518:
	.string	"convertfrom"
.LASF475:
	.string	"mult"
.LASF727:
	.string	"bcol"
.LASF652:
	.string	"lens"
.LASF318:
	.string	"nsends"
.LASF385:
	.string	"MAT_SYMMETRIC"
.LASF86:
	.string	"description"
.LASF462:
	.string	"htype"
.LASF285:
	.string	"mtdot_local"
.LASF608:
	.string	"keepnonzeropattern"
.LASF458:
	.string	"fctx"
.LASF73:
	.string	"PETSC_USE_POINTER"
.LASF600:
	.string	"roworiented"
.LASF424:
	.string	"MAT_GLOBAL_SUM"
.LASF15:
	.string	"unsigned int"
.LASF555:
	.string	"setvaluesrow"
.LASF242:
	.string	"PetscLayout"
.LASF342:
	.string	"factortype"
.LASF432:
	.string	"pivotinblocks"
.LASF174:
	.string	"active"
.LASF373:
	.string	"MAT_IGNORE_MATRIX"
.LASF17:
	.string	"short int"
.LASF587:
	.string	"local_used"
.LASF298:
	.string	"sqrt"
.LASF40:
	.string	"_vtable_offset"
.LASF519:
	.string	"usescaledform"
.LASF212:
	.string	"indices"
.LASF199:
	.string	"curStage"
.LASF508:
	.string	"setblocksize"
.LASF572:
	.string	"multhermitiantransposeadd"
.LASF736:
	.string	"PetscTrMalloc"
.LASF433:
	.string	"zeropivot"
.LASF172:
	.string	"_n_PetscClassPerfLog"
.LASF288:
	.string	"conjugate"
.LASF348:
	.string	"nullsp"
	.ident	"GCC: (GNU) 4.4.6 20110731 (Red Hat 4.4.6-3)"
	.section	.note.GNU-stack,"",@progbits
