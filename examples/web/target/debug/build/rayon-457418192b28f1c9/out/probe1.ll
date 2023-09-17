; ModuleID = 'probe1.1d2d7b3a-cgu.0'
source_filename = "probe1.1d2d7b3a-cgu.0"
target datalayout = "e-m:w-p270:32:32-p271:32:32-p272:64:64-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-pc-windows-msvc"

%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>" = type { %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>" }
%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>" = type { { i32, i32 }, i64, i8, [7 x i8] }

@alloc_4aead6e2018a46d0df208d5729447de7 = private unnamed_addr constant <{ [27 x i8] }> <{ [27 x i8] c"assertion failed: step != 0" }>, align 1
@alloc_ff4ff60397691834652607214c01e676 = private unnamed_addr constant <{ [89 x i8] }> <{ [89 x i8] c"/rustc/90c541806f23a127002de5b4038be731ba1458ca\\library\\core\\src\\iter\\adapters\\step_by.rs" }>, align 1
@alloc_67d1a421b943b48dfd94e6e3778018a5 = private unnamed_addr constant <{ ptr, [16 x i8] }> <{ ptr @alloc_ff4ff60397691834652607214c01e676, [16 x i8] c"Y\00\00\00\00\00\00\00\15\00\00\00\09\00\00\00" }>, align 8

; core::iter::traits::iterator::Iterator::rev
; Function Attrs: inlinehint uwtable
define void @_ZN4core4iter6traits8iterator8Iterator3rev17h2bcafd451d156438E(ptr sret(%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>") %0, ptr %self) unnamed_addr #0 {
start:
; call core::iter::adapters::rev::Rev<T>::new
  call void @"_ZN4core4iter8adapters3rev12Rev$LT$T$GT$3new17hf30b2349873ad53fE"(ptr sret(%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>") %0, ptr %self)
  ret void
}

; core::iter::traits::iterator::Iterator::step_by
; Function Attrs: inlinehint uwtable
define void @_ZN4core4iter6traits8iterator8Iterator7step_by17h3ae1e90a2afb7496E(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") %0, i32 %self.0, i32 %self.1, i64 %step) unnamed_addr #0 {
start:
; call core::iter::adapters::step_by::StepBy<I>::new
  call void @"_ZN4core4iter8adapters7step_by15StepBy$LT$I$GT$3new17hdf44e1e4924cce15E"(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") %0, i32 %self.0, i32 %self.1, i64 %step)
  ret void
}

; core::iter::adapters::rev::Rev<T>::new
; Function Attrs: uwtable
define void @"_ZN4core4iter8adapters3rev12Rev$LT$T$GT$3new17hf30b2349873ad53fE"(ptr sret(%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>") %0, ptr %iter) unnamed_addr #1 {
start:
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %0, ptr align 8 %iter, i64 24, i1 false)
  ret void
}

; core::iter::adapters::step_by::StepBy<I>::new
; Function Attrs: uwtable
define void @"_ZN4core4iter8adapters7step_by15StepBy$LT$I$GT$3new17hdf44e1e4924cce15E"(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") %0, i32 %iter.0, i32 %iter.1, i64 %step) unnamed_addr #1 personality ptr @__CxxFrameHandler3 {
start:
  %_4 = icmp ne i64 %step, 0
  %_3 = xor i1 %_4, true
  br i1 %_3, label %bb1, label %bb2

bb2:                                              ; preds = %start
  %_7 = sub i64 %step, 1
  %1 = getelementptr inbounds { i32, i32 }, ptr %0, i32 0, i32 0
  store i32 %iter.0, ptr %1, align 8
  %2 = getelementptr inbounds { i32, i32 }, ptr %0, i32 0, i32 1
  store i32 %iter.1, ptr %2, align 4
  %3 = getelementptr inbounds %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>", ptr %0, i32 0, i32 1
  store i64 %_7, ptr %3, align 8
  %4 = getelementptr inbounds %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>", ptr %0, i32 0, i32 2
  store i8 1, ptr %4, align 8
  ret void

bb1:                                              ; preds = %start
; invoke core::panicking::panic
  invoke void @_ZN4core9panicking5panic17hb92e8c073af93b07E(ptr align 1 @alloc_4aead6e2018a46d0df208d5729447de7, i64 27, ptr align 8 @alloc_67d1a421b943b48dfd94e6e3778018a5) #5
          to label %unreachable unwind label %funclet_bb3

bb3:                                              ; preds = %funclet_bb3
  cleanupret from %cleanuppad unwind to caller

funclet_bb3:                                      ; preds = %bb1
  %cleanuppad = cleanuppad within none []
  br label %bb3

unreachable:                                      ; preds = %bb1
  unreachable
}

; probe1::probe
; Function Attrs: uwtable
define void @_ZN6probe15probe17h1762f0cda78aa6a3E() unnamed_addr #1 {
start:
  %_3 = alloca { i32, i32 }, align 4
  %_2 = alloca %"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>", align 8
  %_1 = alloca %"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>", align 8
  store i32 0, ptr %_3, align 4
  %0 = getelementptr inbounds { i32, i32 }, ptr %_3, i32 0, i32 1
  store i32 10, ptr %0, align 4
  %1 = getelementptr inbounds { i32, i32 }, ptr %_3, i32 0, i32 0
  %2 = load i32, ptr %1, align 4, !noundef !1
  %3 = getelementptr inbounds { i32, i32 }, ptr %_3, i32 0, i32 1
  %4 = load i32, ptr %3, align 4, !noundef !1
; call core::iter::traits::iterator::Iterator::step_by
  call void @_ZN4core4iter6traits8iterator8Iterator7step_by17h3ae1e90a2afb7496E(ptr sret(%"core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>") %_2, i32 %2, i32 %4, i64 2)
; call core::iter::traits::iterator::Iterator::rev
  call void @_ZN4core4iter6traits8iterator8Iterator3rev17h2bcafd451d156438E(ptr sret(%"core::iter::adapters::rev::Rev<core::iter::adapters::step_by::StepBy<core::ops::range::Range<i32>>>") %_1, ptr %_2)
  ret void
}

; Function Attrs: nocallback nofree nounwind willreturn memory(argmem: readwrite)
declare void @llvm.memcpy.p0.p0.i64(ptr noalias nocapture writeonly, ptr noalias nocapture readonly, i64, i1 immarg) #2

declare i32 @__CxxFrameHandler3(...) unnamed_addr #3

; core::panicking::panic
; Function Attrs: cold noinline noreturn uwtable
declare void @_ZN4core9panicking5panic17hb92e8c073af93b07E(ptr align 1, i64, ptr align 8) unnamed_addr #4

attributes #0 = { inlinehint uwtable "target-cpu"="x86-64" }
attributes #1 = { uwtable "target-cpu"="x86-64" }
attributes #2 = { nocallback nofree nounwind willreturn memory(argmem: readwrite) }
attributes #3 = { "target-cpu"="x86-64" }
attributes #4 = { cold noinline noreturn uwtable "target-cpu"="x86-64" }
attributes #5 = { noreturn }

!llvm.module.flags = !{!0}

!0 = !{i32 8, !"PIC Level", i32 2}
!1 = !{}
